#include "Scheduler.h"

Scheduler::Scheduler()
{
	m_nBound = 0;
	m_nThread = 48;
	m_nGraphBound = 80;
	m_epsilon = 0.1;
	m_nMFLow = 0;
	m_bStop = false;
	m_dPO = 2;

	m_increase = 0.02;
	m_critical = 0.4;
	m_nRun = 30;
}

int Scheduler::CallSMT(NetList& netlist, int nMF)
{//call SMT for sub-netlist
	vector<Node>& vecNode = netlist.m_vecNode;
	vector<Node>& vecIn = netlist.m_vecIn;
	unsigned int nNumNode = vecNode.size();
	unsigned int nTotNumNode = nNumNode + vecIn.size();
	z3::context c;
	z3::optimize s(c);

	z3::set_param("sat.threads", 2);
	z3::params p(c);
	p.set("timeout", (unsigned int)5 * 60 * 1000);
	s.set(p);
	vector<z3::expr_vector> matT;
	for (int n = 0; n < nNumNode; n++) {
		z3::expr_vector tmp(c);
		for (int t = 0; t < nNumNode + 1; t++)
			tmp.push_back(c.bool_const(("T_" + to_string(n) + "," + to_string(t)).c_str()));
		matT.push_back(tmp);
	}

	vector<z3::expr_vector> matA;
	for (int n = 0; n < nTotNumNode; n++) {
		z3::expr_vector tmp(c);
		//tmp.resize(nNumNode + 1);
		matA.push_back(tmp);
	}

	for (unsigned int n = 0; n < nNumNode; n++)
	{
		int nInCone = 0;
		int nOutCone = 0;
		for (int t = 0; t < nInCone + 1; t++)
			s.add(!matT[n][t]);
		for (unsigned int t = nNumNode; t > (nNumNode - 1 - nOutCone); t--)
			s.add(matT[n][t]);
		for (unsigned int t = nInCone + 1; t < (nNumNode - nOutCone - 1); t++)
			s.add(!matT[n][t] || matT[n][t + 1]);
	}

	for (unsigned t = 0; t < nNumNode; t++)
	{
		for (unsigned int i = 1; i < nNumNode; i++)
		{
			for (unsigned int j = 0; j < i; j++)
				s.add(matT[i][t] || !matT[i][t + 1] || matT[j][t] || !matT[j][t + 1]);
		}
	}

	for (Node& Nodei : vecNode)
	{
		unsigned int i = Nodei.m_nIndex;
		for (int j : Nodei.m_vecnSucc)
		{
			for (unsigned t = 0; t < nNumNode; t++)
				s.add(matT[i][t] || !matT[j][t + 1]);
		}
		if (Nodei.m_bPO)
		{
			for (unsigned int t = 0; t < nNumNode + 1; t++)
				matA[i].push_back(matT[i][t]);
		}
		else
		{
			for (unsigned int t = 0; t < nNumNode + 1; t++)
			{
				z3::expr_vector temp(c);
				for (int j : vecNode[i].m_vecnSucc)
					temp.push_back(!matT[j][t]);
				z3::expr tp = mk_or(temp);
				matA[i].push_back((tp && matT[i][t]));
			}
		}
	}

	for (Node& Nodei : vecIn)
	{
		unsigned int i = Nodei.m_nIndex;
		if (Nodei.m_bPO)
		{
			for (unsigned int t = 0; t < nNumNode + 1; t++)
				matA[i].push_back(c.bool_val(true));
		}
		else
		{
			for (unsigned int t = 0; t < nNumNode + 1; t++)
			{
				z3::expr_vector temp(c);
				for (int j : Nodei.m_vecnSucc)
					temp.push_back(!matT[j][t]);
				matA[i].push_back(mk_or(temp));
			}
		}
	}
	z3::expr maxR = c.int_val(0);
	z3::expr_vector vecActInT(c);
	for (unsigned int t = 0; t < nNumNode; t++)
	{
		z3::expr_vector temp(c);
		for (unsigned int n = 0; n < nTotNumNode; n++)
			temp.push_back(matA[n][t + 1]);
		vecActInT.push_back(sum(temp));
	}

	if (nMF == 0)
	{
		for (unsigned int t = 0; t < nNumNode; t++)
			maxR = max(vecActInT[t], maxR);
		s.minimize(maxR);

	}
	else
	{
		for (unsigned int t = 0; t < nNumNode; t++)
			s.add(vecActInT[t] <= nMF);
	}
	int nRet = -1;
	//m_netlist.m_vecnSchedule.clear();
	if (netlist.m_vecnSchedule.empty())
		netlist.m_vecnSchedule.assign(nNumNode, 0);
	/*
	if (nMF == 0)
		cout << "Running minimize mode\n";
	else
		cout << "Runnning feasibility mode, " << "Checking for " << nMF << " rows\n";

	clock_t start = clock();
	if (s.check() != z3::sat)
		cout << "Infeasible for " << nMF << " rows\n";
	else
	*/
	if (s.check() == z3::sat)
	{
		z3::model m = s.get_model();
		for (int n = 0; n < nNumNode; n++)
		{
			int nSched = nNumNode - m.eval(sum(matT[n])).get_numeral_int();
			//cout << "S" << n << ": " << nSched << "\n";
			netlist.m_vecnSchedule[nSched] = n;
		}
		for (unsigned int t = 0; t < nNumNode; t++)
			nRet = max(m.eval(vecActInT[t]).get_numeral_int(), nRet);
		//cout << "Feasible for " << nRet << " rows\n";
	}
	//cout << "Use: " << (clock() - start) * 1000 / double(CLOCKS_PER_SEC) << " ms\n";
	return nRet;
}

int Scheduler::BiSMT(NetList& netlist, int nUpper, int nLower)
{//binary search with SMT for sub-netlist
	nLower = netlist.m_vecnPO.size();//tmp
	int nRes = CallSMT(netlist, nLower);
	if (nRes != -1)
		return nRes;

	nRes = CallSMT(netlist, nUpper);
	if (nRes == -1)
		return -1;
	nUpper = nRes;

	while ((nUpper - nLower) > 1)
	{
		if (m_bStop)
			break;
		if (m_nMFLow >= nUpper)//tmp
			break;

		int nPres = round((nUpper + nLower) / 2.0);
		nRes = CallSMT(netlist, nPres);
		if (nRes == -1)
			nLower = nPres;
	//nLower = max(nPres, m_nMFLow - 1);//tmp
		else
			nUpper = nRes;
	}
	return nUpper;
}

int Scheduler::PartitionILP(NetList& netlist, double dPOWeight)
{//ILP-based partition
	vector<Node>& vecNode = netlist.m_vecNode;
	vector<Node>& vecIn = netlist.m_vecIn;
	unsigned int nNumNode = (unsigned int)vecNode.size();
	cout << "Partitioning " << nNumNode << " nodes with weight = " << dPOWeight << "\n";;

	unsigned int nIn = (unsigned int)vecIn.size();
	GRBEnv env = GRBEnv(true);
	env.set("LogFile", "mip1.log");
	env.start();
	GRBModel model = GRBModel(env);
	model.set(GRB_IntParam_OutputFlag, 0);
	//model.set(GRB_DoubleParam_MIPGap, 0.1);
	model.set(GRB_DoubleParam_TimeLimit, 3 * 60);
	//model.set(GRB_IntParam_Threads, m_nThread);
	GRBVar* vecP = new GRBVar[nNumNode];
	GRBVar* vecSucc = new GRBVar[nNumNode + nIn];
	for (unsigned int n = 0; n < nNumNode; n++)
	{
		vecP[n] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "P_" + to_string(n));
		vecSucc[n] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "Succ_" + to_string(n));
		//if (n >= nNumNode)
			//model.addConstr(vecP[n] == 0);
	}
	for (unsigned int n = nNumNode; n < nNumNode + nIn; n++)
	{
		vecSucc[n] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "Succ_" + to_string(n));
	}

	for (Node& Nodei : vecNode)
	{
		unsigned int i = Nodei.m_nIndex;
		vector<GRBVar> temp;
		for (int j : Nodei.m_vecnSucc)
		{
			model.addConstr(vecP[i] <= vecP[j]);
			temp.push_back(vecP[j]);
		}
		if (Nodei.m_bPO)
			model.addConstr(vecSucc[i] == 1);
		else
			model.addGenConstrOr(vecSucc[i], temp.data(), (int)temp.size());
	}
	for (Node& Nodei : vecIn)
	{
		unsigned int i = Nodei.m_nIndex;
		if (Nodei.m_bPO)
			model.addConstr(vecSucc[i] == 1);
		else
		{
			vector<GRBVar> temp;
			for (int j : Nodei.m_vecnSucc)
				temp.push_back(vecP[j]);
			model.addGenConstrOr(vecSucc[i], temp.data(), (int)temp.size());
		}
	}

	GRBLinExpr expr = 0;
	for (unsigned int n = 0; n < nNumNode; n++)
		expr += vecP[n];
	model.addConstr(expr, GRB_LESS_EQUAL, (1 + m_epsilon) * (double)nNumNode / 2.0);
	model.addConstr(expr, GRB_GREATER_EQUAL, (1 - m_epsilon) * (double)nNumNode / 2.0);
	GRBQuadExpr qexpr = 0;
	for (unsigned int n = 0; n < nNumNode; n++)
	{
		//if (vecNode[n].m_bPO && (vecNode[n].m_vecnSucc.size() == 0))
		if (vecNode[n].m_bPO)
			qexpr += (1 - vecP[n]) * vecSucc[n] * std::max(1.0, dPOWeight);
		else
			qexpr += (1 - vecP[n]) * vecSucc[n];
	}
	for (unsigned int n = nNumNode; n < nNumNode + nIn; n++)
		qexpr += vecSucc[n];

	model.setObjective(qexpr, GRB_MINIMIZE);

	int nRet = -1;
	int nSecondHalf = 0;
	int nFront = 0;
	clock_t start = clock();
	model.optimize();
	if (model.get(GRB_IntAttr_Status) == GRB_INFEASIBLE)
		cout << "Partition BUG\n";
	else
	{
		for (unsigned int n = 0; n < nNumNode; n++)
		{
			//cout << vecP[n].get(GRB_StringAttr_VarName) << " "
				//<< (unsigned int)round(vecP[n].get(GRB_DoubleAttr_X)) << "\n";
			if (round(vecP[n].get(GRB_DoubleAttr_X)) == 1)
			{
				vecNode[n].m_nPartition = 1;
				nSecondHalf++;
			}
			else if (round(vecSucc[n].get(GRB_DoubleAttr_X)) == 0)
				vecNode[n].m_nPartition = 0;
			else
			{
				vecNode[n].m_nPartition = 2;
				nFront++;
			}
		}
		for (unsigned int n = 0; n < nIn; n++)
		{
			if (round(vecSucc[n + nNumNode].get(GRB_DoubleAttr_X)) == 0)
				vecIn[n].m_nPartition = 0;
			else
			{
				vecIn[n].m_nPartition = 2;
				nFront++;
			}
		}
		cout << "Partitioned into: " << nNumNode - nSecondHalf << " : " << nSecondHalf << "\n";
		cout << "#Frontier=" << nFront << "\n";
		nRet = nFront;
	}
	cout << "Use: " << double(clock() - start) * 1000 / CLOCKS_PER_SEC << " ms\n";

	delete[] vecP;
	delete[] vecSucc;
	return nRet;
}

void Scheduler::IterPart(NetList& netlist, double dLeft, double dRight)
{//iterative partition
	vector<Node>& vecNode = netlist.m_vecNode;
	vector<Node>& vecIn = netlist.m_vecIn;
	vector<int>& vecnPo = netlist.m_vecnPO;
	if (vecnPo.size() > m_nMFLow)
		m_nMFLow = vecnPo.size();
	if (vecnPo.size() > m_nBound)
		m_bStop = true;
	if (m_bStop)
		return;
	int nTotSize = vecNode.size();
	if (nTotSize <= m_nGraphBound)
	{
		m_vecPartNet.push_back(netlist);
		return;
	}
	//double dPOWeight = sqrt(dLeft * dRight);
	double dPOWeight = ((double)dLeft - 1.0) / 2 + 1;
	//double dPOWeight = sqrt(dLeft);
	//double dPOWeight = std::max(1.0, 0.9 * dLeft);
	PartitionILP(netlist, dLeft);

	NetList netlist0, netlist1;
	vector<Node>& vecNode0 = netlist0.m_vecNode;
	vector<Node>& vecNode1 = netlist1.m_vecNode;
	vector<Node>& vecIn0 = netlist0.m_vecIn;
	vector<Node>& vecIn1 = netlist1.m_vecIn;
	vector<int>& vecnNewPo0 = netlist0.m_vecnPO;
	vector<int>& vecnNewPo1 = netlist1.m_vecnPO;

	vector<int> vnPart0;
	vector<int> vnPart1;
	vector<int> vnFront;
	for (unsigned int n = 0; n < vecNode.size(); n++)
	{
		int nPart = vecNode[n].m_nPartition;
		if (nPart == 0)
			vnPart0.push_back(n);
		else if (nPart == 1)
			vnPart1.push_back(n);
		else
		{
			vnPart0.push_back(n);
			vnFront.push_back(n);
		}
	}
	sort(vnPart0.begin(), vnPart0.end());
	sort(vnPart1.begin(), vnPart1.end());
	sort(vnFront.begin(), vnFront.end());
	std::map<int, int> map0OldNew;
	std::map<int, int> map1OldNew;
	std::map<int, int> mapFOldNew;
	int nPart0Size = vnPart0.size();
	int nPart1Size = vnPart1.size();
	int nFrontier = vnFront.size();
	for (unsigned int n = 0; n < nPart0Size; n++)
		map0OldNew[vnPart0[n]] = n;
	for (unsigned int n = 0; n < nPart1Size; n++)
		map1OldNew[vnPart1[n]] = n;
	for (unsigned int n = 0; n < nFrontier; n++)
		mapFOldNew[vnFront[n]] = n;

	vector<Node> tmp0(nPart0Size);
	vecNode0 = tmp0;
	for (unsigned int i = 0; i < nPart0Size; i++)
	{
		vecNode0[i].m_nIndex = i;
		vecNode0[i].m_nOrigIndex = vecNode[vnPart0[i]].m_nOrigIndex;
		for (int jold : vecNode[vnPart0[i]].m_vecnSucc)
		{
			if (map0OldNew.find(jold) != map0OldNew.end())
				vecNode0[i].m_vecnSucc.push_back(map0OldNew[jold]);
		}
	}

	vecnNewPo0.clear();
	for (int front : vnFront)
	{
		vecnNewPo0.push_back(map0OldNew[front]);
		vecNode0[map0OldNew[front]].m_bPO = true;
	}

	vector<Node> tmp0i(vecIn.size());
	vecIn0 = tmp0i;
	for (unsigned int i = 0; i < vecIn.size(); i++)
	{
		vecIn0[i].m_nIndex = nPart0Size + i;
		vecIn0[i].m_nOrigIndex = vecIn[i].m_nOrigIndex;
		for (int jold : vecIn[i].m_vecnSucc)
		{
			if (map0OldNew.find(jold) != map0OldNew.end())
				vecIn0[i].m_vecnSucc.push_back(map0OldNew[jold]);
		}
		if (vecIn[i].m_nPartition == 2)
		{
			vecnNewPo0.push_back(nPart0Size + i);
			vecIn0[i].m_bPO = true;
		}
	}
	//IterPart(netlist0, dLeft, dPOWeight);
	IterPart(netlist0, dPOWeight, dRight);

	vector<Node> tmp1(nPart1Size);
	vecNode1 = tmp1;
	for (unsigned int i = 0; i < nPart1Size; i++)
	{
		vecNode1[i].m_nIndex = i;
		vecNode1[i].m_nOrigIndex = vecNode[vnPart1[i]].m_nOrigIndex;
		for (int jold : vecNode[vnPart1[i]].m_vecnSucc)
		{
			if (map1OldNew.find(jold) != map1OldNew.end())
				vecNode1[i].m_vecnSucc.push_back(map1OldNew[jold]);
		}
	}
	vector<int> vecnInFront;
	for (unsigned int i = 0; i < vecIn.size(); i++)
	{
		if (vecIn[i].m_nPartition == 2)
			vecnInFront.push_back(vecIn[i].m_nIndex);
	}
	int nInFront = vecnInFront.size();

	vector<Node> tmp1i(nFrontier + nInFront);
	vecIn1 = tmp1i;
	for (unsigned int i = 0; i < nFrontier; i++)
	{
		vecIn1[i].m_nIndex = nPart1Size + i;
		vecIn1[i].m_nOrigIndex = vecNode[vnFront[i]].m_nOrigIndex;
		for (int jold : vecNode[vnFront[i]].m_vecnSucc)
		{
			if (map1OldNew.find(jold) != map1OldNew.end())
				vecIn1[i].m_vecnSucc.push_back(map1OldNew[jold]);
		}
	}

	vecnNewPo1.clear();
	for (int po : vecnPo)
	{
		if (map1OldNew.find(po) != map1OldNew.end())
		{
			vecnNewPo1.push_back(map1OldNew[po]);
			vecNode1[map1OldNew[po]].m_bPO = true;
		}
		else if (mapFOldNew.find(po) != mapFOldNew.end())
		{
			vecnNewPo1.push_back(mapFOldNew[po] + nPart1Size);
			vecIn1[mapFOldNew[po]].m_bPO = true;
		}
	}
	for (unsigned int i = 0; i < nInFront; i++)
	{
		vecIn1[i + nFrontier].m_nIndex = nPart1Size + nFrontier + i;
		vecIn1[i + nFrontier].m_nOrigIndex = vecIn[vecnInFront[i] - nTotSize].m_nOrigIndex;
		for (int jold : vecIn[vecnInFront[i] - nTotSize].m_vecnSucc)
		{
			if (map1OldNew.find(jold) != map1OldNew.end())
				vecIn1[i + nFrontier].m_vecnSucc.push_back(map1OldNew[jold]);
		}
		if (vecIn[vecnInFront[i] - nTotSize].m_bPO)
		{
			vecnNewPo1.push_back(nPart1Size + nFrontier + i);
			vecIn1[i + nFrontier].m_bPO = true;
		}
	}
	IterPart(netlist1, dPOWeight, dRight);
}

void Scheduler::ScheduleThread()
{//scheduling thread
	while (WaitForSingleObject(m_hEventKillThread, 1) != WAIT_OBJECT_0)
	{
		if (m_bStop)
			break;
		LONG nID = InterlockedIncrement(&m_nProcNet);
		if (nID > m_nTotNet)
			break;
		int nSch = BiSMT(m_vecPartNet[nID - 1], m_nBound, m_nMFLow);
		if (nSch == -1)
		{
			m_bStop = true;
			break;
		}
		EnterCriticalSection(&m_rCritical);
		if (m_nMFLow < nSch)
			m_nMFLow = nSch;
		m_nProgressDone += m_vecPartNet[nID - 1].m_vecNode.size();
		LeaveCriticalSection(&m_rCritical);
		cout << "=====================Progress: scheduled " << m_nProgressDone << "/" << m_nProgressTotal << " nodes. Sub MF: "
			<< nSch << "\n";
	}
}

DWORD WINAPI MyThreadFunction(LPVOID lpParam)
{
	if (lpParam == NULL)
		return 0;

	Scheduler* pMyScheduler = (Scheduler*)lpParam;
	pMyScheduler->ScheduleThread();
	return 0;
}

void Scheduler::ThreadIterPartScheduler()
{//multi-thread interative partition scheduler
	cout << "Scheduling netlist with " << m_netlist.m_vecNode.size() << " nodes with MF bound " << m_nBound << "\n";
	m_bStop = false;
	m_nProgressTotal = m_netlist.m_vecNode.size();
	m_nProgressDone = 0;
	m_nMFLow = 0;
	m_vecPartNet.clear();
	clock_t start = clock();
	IterPart(m_netlist, m_dPO, 1);
	if (m_bStop)
	{
		cout << "MF after partition too large\n";
		return;
	}
	m_nTotNet = m_vecPartNet.size();
	m_nProcNet = 0;
	cout << "Partition Run time: " << double(clock() - start) * 1000 / CLOCKS_PER_SEC << " ms\n";
	cout << "MF lower bound: " << m_nMFLow << "\n";

	InitializeCriticalSection(&m_rCritical);
	m_hEventKillThread = CreateEvent(NULL, TRUE, FALSE, NULL); // manual reset, initially reset
	HANDLE* pThreadHandle = new HANDLE[m_nThread];
	PDWORD pdwThreadID = new DWORD[m_nThread];
	for (int i = 0; i < m_nThread; i++)
	{
		pThreadHandle[i] = CreateThread(NULL, 0, MyThreadFunction, this, 0, pdwThreadID + i);
		assert(pThreadHandle[i] != NULL);
	}
	WaitForMultipleObjects(m_nThread, pThreadHandle, TRUE, INFINITE);
	for (int i = 0; i < m_nThread; i++)
		CloseHandle(pThreadHandle[i]);

	delete[] pThreadHandle;
	delete[] pdwThreadID;
	CloseHandle(m_hEventKillThread);
	DeleteCriticalSection(&m_rCritical);

	for (NetList& net : m_vecPartNet)
	{
		for (int n : net.m_vecnSchedule)
			m_netlist.m_vecnSchedule.push_back(net.m_vecNode[n].m_nOrigIndex);
	}
	m_nFootPrint = m_nMFLow;
	if (m_nFootPrint < m_netlist.m_vecIn.size())
		m_nFootPrint = m_netlist.m_vecIn.size();
	cout << "Schedule ends, MF = " << m_nFootPrint << "\n";
}

bool Scheduler::IsDivisor(const Node& ndRoot, const Node& ndDiv)
{//whether ndDiv can be divisor for ndRoot
	const set<int>& setRootConePI = ndRoot.m_setnConePI;
	int nNotIn = 0;
	for (int pi : ndDiv.m_setnConePI)
	{
		if (setRootConePI.find(pi) == setRootConePI.end())
			nNotIn++;
	}
	if (nNotIn < 2)
		return true;
	return false;
}

void Scheduler::ResubThread()
{//thread for MF-Resub
	while (WaitForSingleObject(m_hEventKillThread, 1) != WAIT_OBJECT_0)
	{
		if (m_bFoundResub)
			break;
		LONG m = InterlockedIncrement(&m_nCheckPeak);
		if (m >= m_vecPeakIndex.size())
			break;
		clock_t start = clock();
		NetList netlist = m_curNetlist;
		netlist.m_net = m_curNetlist.m_net.clone();
		int nPeakEnd = m_vecPeakIndex.back();
		Node& nd = netlist.m_vecNode[m_vecPeakIndex[m]];
		if (nd.m_bPO)
			continue;
		int nFanOut = -1;
		for (int fo : nd.m_vecnSucc)
		{
			if (fo > nPeakEnd)
			{
				if (nFanOut != -1)
				{
					nFanOut = -1;
					break;
				}
				nFanOut = fo;
			}
		}
		if (nFanOut == -1)
			continue;
		int foID = nFanOut + netlist.m_nOffset;
		int mID = m_vecPeakIndex[m] + netlist.m_nOffset;
		//cout << "node " << mID << " has single FO " << foID << "\n";
		vector<int> vecCand;
		for (int p : m_vecPeakIndex)
		{
			if (p != m_vecPeakIndex[m] && IsDivisor(netlist.m_vecNode[nFanOut], netlist.m_vecNode[p]))
				vecCand.push_back(p + netlist.m_nOffset);
		}
		for (int i = nPeakEnd + 1; i < nFanOut; i++)
		{
			if (IsDivisor(netlist.m_vecNode[nFanOut], netlist.m_vecNode[i]))
				vecCand.push_back(i + netlist.m_nOffset);
		}

		for (int pi : netlist.m_vecNode[nFanOut].m_setnConePI)
			vecCand.push_back(pi + 1);
		vecCand.push_back(0);

		int nCandSize = vecCand.size();
		//cout << "#candidate nodes: " << nCandSize << "\n";
		if (nCandSize < 3)
			continue;
		InterlockedIncrement(&m_nRoot);
		auto tt = m_tts[foID];
		bool bRand = true;
		//int nTrial = 150000;
		//if (nCandSize <= 100)
		int nTrial = 1500000;
		if (nCandSize <= 200)
		{
			bRand = false;
			nTrial = nCandSize * nCandSize * nCandSize;
		}
		for (int t = 0; t < nTrial; t++)
		{
			if (m_bFoundResub)
				break;
			int xID, yID, zID;
			if (bRand)
			{
				xID = vecCand[rand() % nCandSize];
				while (true)
				{
					yID = vecCand[rand() % nCandSize];
					if (yID != xID)
						break;
				}
				while (true)
				{
					zID = vecCand[rand() % nCandSize];
					if (zID != xID && zID != yID)
						break;
				}

			}
			else
			{
				xID = vecCand[t % nCandSize];
				yID = vecCand[(t / nCandSize) % nCandSize];
				zID = vecCand[t / nCandSize / nCandSize];
				if (xID <= yID || yID <= zID)
					continue;
			}
			validator_params vps;
			vps.max_clauses = 1000;
			vps.conflict_limit = 100;

			const kitty::partial_truth_table& ttx = m_tts[xID];
			const kitty::partial_truth_table& tty = m_tts[yID];
			const kitty::partial_truth_table& ttz = m_tts[zID];

			vector<int> vecnPredID = { xID, yID, zID };
			vector<bool> vecbPredComp = { m_vecPhase[xID] , m_vecPhase[yID] , m_vecPhase[zID] };
			if (tt == detail::ternary_xor(ttx, tty, ttz))
			{
				xmg_network ntk_s = netlist.m_net.clone();
				xmg_network::signal sx = m_vecPhase[xID] ? !ntk_s.make_signal(xID) : ntk_s.make_signal(xID);
				xmg_network::signal sy = m_vecPhase[yID] ? !ntk_s.make_signal(yID) : ntk_s.make_signal(yID);
				xmg_network::signal sz = m_vecPhase[zID] ? !ntk_s.make_signal(zID) : ntk_s.make_signal(zID);
				circuit_validator<xmg_network, bill::solvers::bsat2> validator_s(ntk_s, vps);
				//cout << "--------------------------------------Find candidate XOR " << xID << " " << yID << " " << zID << "\n";
				xmg_network::signal candidate = m_vecPhase[foID] ? !ntk_s.create_xor3(sx, sy, sz) : ntk_s.create_xor3(sx, sy, sz);
				optional<bool> res = validator_s.validate(ntk_s.make_signal(foID), candidate);
				if (!res || !(*res))
				{
					InterlockedIncrement(&m_nVFail);
					continue;
				}
				else
				{
					EnterCriticalSection(&m_rCritical);
					//cout << "----------Substitute fanins of " << foID << " into XOR " << xID << " " << yID << " " << zID << "\n";
					//cout << vecbPredComp[0] << " " << vecbPredComp[1] << " " << vecbPredComp[2] << " " << m_vecPhase[foID] << "\n";
					netlist.SubPred(nFanOut, vecnPredID, vecbPredComp, false, m_vecPhase[foID]);
					m_newNetlist = netlist;
					m_bFoundResub = true;
					LeaveCriticalSection(&m_rCritical);
					//ntk_s.substitute_node(foID, candidate);
					//write_verilog(ntk_s, "ttt.v");
					//netlist.m_net = ntk_s.clone();
					//cout << "Valid Fail: " << m_nVFail << "\n";
					break;
				}
			}
			if (tt == kitty::ternary_majority(ttx, tty, ttz))
			{
				xmg_network ntk_s = netlist.m_net.clone();
				xmg_network::signal sx = m_vecPhase[xID] ? !ntk_s.make_signal(xID) : ntk_s.make_signal(xID);
				xmg_network::signal sy = m_vecPhase[yID] ? !ntk_s.make_signal(yID) : ntk_s.make_signal(yID);
				xmg_network::signal sz = m_vecPhase[zID] ? !ntk_s.make_signal(zID) : ntk_s.make_signal(zID);
				circuit_validator<xmg_network, bill::solvers::bsat2> validator_s(ntk_s, vps);
				//cout << "--------------------------------------Find candidate MAJ " << xID << " " << yID << " " << zID << "\n";
				xmg_network::signal candidate = m_vecPhase[foID] ? !ntk_s.create_maj(sx, sy, sz) : ntk_s.create_maj(sx, sy, sz);
				optional<bool> res = validator_s.validate(ntk_s.make_signal(foID), candidate);
				if (!res || !(*res))
				{
					InterlockedIncrement(&m_nVFail);
					continue;
				}
				else
				{
					EnterCriticalSection(&m_rCritical);
					//cout << "----------Substitute fanins of " << foID << " into MAJ " << xID << " " << yID << " " << zID << "\n";
					//cout << vecbPredComp[0] << " " << vecbPredComp[1] << " " << vecbPredComp[2] << " " << m_vecPhase[foID] << "\n";
					netlist.SubPred(nFanOut, vecnPredID, vecbPredComp, true, m_vecPhase[foID]);
					m_newNetlist = netlist;
					m_bFoundResub = true;
					LeaveCriticalSection(&m_rCritical);
					//ntk_s.substitute_node(foID, candidate);
					//write_verilog(ntk_s, "ttt.v");
					//netlist.m_net = ntk_s.clone();
					//cout << "Valid Fail: " << m_nVFail << "\n";
					break;
				}
			}
			if (tt == kitty::ternary_majority(~ttx, tty, ttz))
			{
				xmg_network ntk_s = netlist.m_net.clone();
				xmg_network::signal sx = m_vecPhase[xID] ? !ntk_s.make_signal(xID) : ntk_s.make_signal(xID);
				xmg_network::signal sy = m_vecPhase[yID] ? !ntk_s.make_signal(yID) : ntk_s.make_signal(yID);
				xmg_network::signal sz = m_vecPhase[zID] ? !ntk_s.make_signal(zID) : ntk_s.make_signal(zID);
				circuit_validator<xmg_network, bill::solvers::bsat2> validator_s(ntk_s, vps);
				//cout << "--------------------------------------Find candidate NMAJ " << xID << " " << yID << " " << zID << "\n";
				xmg_network::signal candidate = m_vecPhase[foID] ? !ntk_s.create_maj(!sx, sy, sz) : ntk_s.create_maj(!sx, sy, sz);
				optional<bool> res = validator_s.validate(ntk_s.make_signal(foID), candidate);
				if (!res || !(*res))
				{
					InterlockedIncrement(&m_nVFail);
					continue;
				}
				else
				{
					EnterCriticalSection(&m_rCritical);
					//cout << "----------Substitute fanins of " << foID << " into NMAJ0 " << xID << " " << yID << " " << zID << "\n";
					vecbPredComp[0] = !vecbPredComp[0];
					//cout << vecbPredComp[0] << " " << vecbPredComp[1] << " " << vecbPredComp[2] << " " << m_vecPhase[foID] << "\n";
					netlist.SubPred(nFanOut, vecnPredID, vecbPredComp, true, m_vecPhase[foID]);
					m_newNetlist = netlist;
					m_bFoundResub = true;
					LeaveCriticalSection(&m_rCritical);
					//ntk_s.substitute_node(foID, candidate);
					//write_verilog(ntk_s, "ttt.v");
					//netlist.m_net = ntk_s.clone();
					//cout << "Valid Fail: " << m_nVFail << "\n";
					break;
				}
			}
			if (tt == kitty::ternary_majority(ttx, ~tty, ttz))
			{
				xmg_network ntk_s = netlist.m_net.clone();
				xmg_network::signal sx = m_vecPhase[xID] ? !ntk_s.make_signal(xID) : ntk_s.make_signal(xID);
				xmg_network::signal sy = m_vecPhase[yID] ? !ntk_s.make_signal(yID) : ntk_s.make_signal(yID);
				xmg_network::signal sz = m_vecPhase[zID] ? !ntk_s.make_signal(zID) : ntk_s.make_signal(zID);
				circuit_validator<xmg_network, bill::solvers::bsat2> validator_s(ntk_s, vps);
				//cout << "--------------------------------------Find candidate NMAJ " << xID << " " << yID << " " << zID << "\n";
				xmg_network::signal candidate = m_vecPhase[foID] ? !ntk_s.create_maj(sx, !sy, sz) : ntk_s.create_maj(sx, !sy, sz);
				optional<bool> res = validator_s.validate(ntk_s.make_signal(foID), candidate);
				if (!res || !(*res))
				{
					InterlockedIncrement(&m_nVFail);
					continue;
				}
				else
				{
					EnterCriticalSection(&m_rCritical);
					//cout << "----------Substitute fanins of " << foID << " into NMAJ1 " << xID << " " << yID << " " << zID << "\n";
					vecbPredComp[1] = !vecbPredComp[1];
					//cout << vecbPredComp[0] << " " << vecbPredComp[1] << " " << vecbPredComp[2] << " " << m_vecPhase[foID] << "\n";
					netlist.SubPred(nFanOut, vecnPredID, vecbPredComp, true, m_vecPhase[foID]);
					m_newNetlist = netlist;
					m_bFoundResub = true;
					LeaveCriticalSection(&m_rCritical);
					//ntk_s.substitute_node(foID, candidate);
					//write_verilog(ntk_s, "ttt.v");
					//netlist.m_net = ntk_s.clone();
					//cout << "Valid Fail: " << m_nVFail << "\n";
					break;
				}
			}
			if (tt == kitty::ternary_majority(ttx, tty, ~ttz))
			{
				xmg_network ntk_s = netlist.m_net.clone();
				xmg_network::signal sx = m_vecPhase[xID] ? !ntk_s.make_signal(xID) : ntk_s.make_signal(xID);
				xmg_network::signal sy = m_vecPhase[yID] ? !ntk_s.make_signal(yID) : ntk_s.make_signal(yID);
				xmg_network::signal sz = m_vecPhase[zID] ? !ntk_s.make_signal(zID) : ntk_s.make_signal(zID);
				circuit_validator<xmg_network, bill::solvers::bsat2> validator_s(ntk_s, vps);
				//cout << "--------------------------------------Find candidate NMAJ " << xID << " " << yID << " " << zID << "\n";
				xmg_network::signal candidate = m_vecPhase[foID] ? !ntk_s.create_maj(sx, sy, !sz) : ntk_s.create_maj(sx, sy, !sz);
				optional<bool> res = validator_s.validate(ntk_s.make_signal(foID), candidate);
				if (!res || !(*res))
				{
					InterlockedIncrement(&m_nVFail);
					continue;
				}
				else
				{
					EnterCriticalSection(&m_rCritical);
					//cout << "----------Substitute fanins of " << foID << " into NMAJ2 " << xID << " " << yID << " " << zID << "\n";
					vecbPredComp[2] = !vecbPredComp[2];
					//cout << vecbPredComp[0] << " " << vecbPredComp[1] << " " << vecbPredComp[2] << " " << m_vecPhase[foID] << "\n";
					netlist.SubPred(nFanOut, vecnPredID, vecbPredComp, true, m_vecPhase[foID]);
					m_newNetlist = netlist;
					m_bFoundResub = true;
					LeaveCriticalSection(&m_rCritical);
					//ntk_s.substitute_node(foID, candidate);
					//write_verilog(ntk_s, "ttt.v");
					//netlist.m_net = ntk_s.clone();
					//cout << "Valid Fail: " << m_nVFail << "\n";
					break;
				}
			}
		}
		//cout << mID<< " #candidate nodes: " << nCandSize << " thread time: " << double(clock() - start) / CLOCKS_PER_SEC << " s\n";
	}
}


DWORD WINAPI ThreadResub(LPVOID lpParam)
{
	if (lpParam == NULL)
		return 0;

	Scheduler* pMyScheduler = (Scheduler*)lpParam;
	pMyScheduler->ResubThread();
	return 0;
}


bool Scheduler::ResubMult(NetList& netlist)
{//multi-thread resub
	m_nVFail = 0;
	m_nRoot = 0;
	xmg_network ntk = netlist.m_net.clone();
	int nNumNode = netlist.m_vecNode.size();
	m_vecPeakIndex = netlist.GetPeak();
	int nPeakSize = m_vecPeakIndex.size();
	//cout << "peak with " << nPeakSize << " nodes\n";

	//simulate
	partial_simulator sim(ntk.num_pis(), 1000);
	pattern_generation_params ps;
	ps.odc_levels = 5;
	pattern_generation(ntk, sim, ps);
	unordered_node_map<kitty::partial_truth_table, xmg_network> tts(ntk);
	simulate_nodes<xmg_network>(ntk, tts, sim, true);
	validator_params vps;
	vps.max_clauses = 1000;
	vps.conflict_limit = 100;
	circuit_validator<xmg_network, bill::solvers::bsat2> validator(ntk, vps);
	vector<bool> vecPhase(tts.size());
	for (int n = 0; n < tts.size(); n++)
	{
		auto& tt = tts[n];
		if (kitty::get_bit(tt, 0))
		{
			tt = ~tt;
			vecPhase[n] = true;
		}
		else
			vecPhase[n] = false;
	}

	//cout << "Checking for resub in peak\n";
	for (int m = 0; m < nPeakSize; m++)
	{
		int mID = m_vecPeakIndex[m] + netlist.m_nOffset;
		/*
		if (tts[mID] == tts[0])
		{
			cout << "node " << mID << "candidate constant\n";
			xmg_network::signal candidate = vecPhase[mID] ? ntk.get_constant(true) : ntk.get_constant(false);
			optional<bool> res = validator.validate(ntk.make_signal(mID), candidate);
			if (!res || !(*res))
				cout << "checked false \n";
			else
			{
				cout << "checked true!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
				ntk.substitute_node(mID, candidate);
				return true;
			}
		}
		*/
		for (int n = m + 1; n < nPeakSize; n++)
		{
			int nID = m_vecPeakIndex[n] + netlist.m_nOffset;
			if (tts[mID] == tts[nID])
			{
				//cout << "node " << mID << " candidate equivalent to node " << nID << "\n";
				xmg_network::signal candidate = (vecPhase[mID] ^ vecPhase[nID]) ? !ntk.make_signal(nID) : ntk.make_signal(nID);
				optional<bool> res = validator.validate(ntk.make_signal(mID), candidate);
				if (!res || !(*res))
					continue;
				else
				{
					//cout << "----------Substitute node " << mID << " with " << nID << " " << (vecPhase[mID] ^ vecPhase[nID]) << "\n";
					netlist.SubNode(std::min(m_vecPeakIndex[m], m_vecPeakIndex[n]), std::max(m_vecPeakIndex[m], m_vecPeakIndex[n]), vecPhase[mID] ^ vecPhase[nID]);
					return true;
				}
			}
		}
	}

	//cout << "Checking for resub in fanout of peak\n";
	m_bFoundResub = false;
	m_curNetlist = netlist;
	m_nCheckPeak = -1;
	m_tts.clear();
	for (int n = 0; n < tts.size(); n++)
		m_tts.push_back(tts[n]);
	m_vecPhase = vecPhase;
	InitializeCriticalSection(&m_rCritical);
	m_hEventKillThread = CreateEvent(NULL, TRUE, FALSE, NULL); // manual reset, initially reset
	HANDLE* pThreadHandle = new HANDLE[m_nThread];
	PDWORD pdwThreadID = new DWORD[m_nThread];
	for (int i = 0; i < m_nThread; i++)
	{
		pThreadHandle[i] = CreateThread(NULL, 0, ThreadResub, this, 0, pdwThreadID + i);
		assert(pThreadHandle[i] != NULL);
	}
	WaitForMultipleObjects(m_nThread, pThreadHandle, TRUE, INFINITE);
	for (int i = 0; i < m_nThread; i++)
		CloseHandle(pThreadHandle[i]);

	delete[] pThreadHandle;
	delete[] pdwThreadID;
	CloseHandle(m_hEventKillThread);
	DeleteCriticalSection(&m_rCritical);

	//cout << "===== root = " << m_nRoot << ", valid fail = " << m_nVFail << "\n";
	if (m_bFoundResub)
	{
		netlist = m_newNetlist;
		return true;
	}
	//cout << "no resub\n";
	return false;
}

bool Scheduler::UpdatePareto(NetList& NewNetList, vector<NetList>& vecParetoList)
{//update pareto
	int nNewSize = NewNetList.m_nSize;
	int nNewMF = NewNetList.m_nMF;
	bool bPareto = true;
	for (auto nb = vecParetoList.begin(); nb != vecParetoList.end();) {
		const NetList& paretopt = (*nb);
		if ((paretopt.m_nSize < nNewSize && paretopt.m_nMF <= nNewMF)
			|| (paretopt.m_nSize <= nNewSize && paretopt.m_nMF < nNewMF))
		{
			bPareto = false;
			break;
		}
		if ((nNewSize < paretopt.m_nSize && nNewMF <= paretopt.m_nMF)
			|| (nNewSize <= paretopt.m_nSize && nNewMF < paretopt.m_nMF))
			nb = vecParetoList.erase(nb);
		else
			++nb;
	}
	if (bPareto)
		vecParetoList.push_back(NewNetList);
	return bPareto;
}

void Scheduler::NewDSE()
{
	clock_t start = clock();
	int nOrigOffset = m_netlist.m_nOffset;
	int nAcceptNode = m_netlist.m_vecNode.size() * (1 + m_increase);
	std::ofstream ofCSV;
	string strResultName = m_netlist.m_strBench + "_result.csv";
	ofCSV.open(strResultName.c_str(), std::ofstream::out);
	cout << "++++++++++++++++++++Bench: " << m_netlist.m_strBench << "++++++++++++++++++++\n";
	cout << "-----------------Round 0:-----------------\n";
	ThreadIterPartScheduler();
	NetList ReorderedNet = m_netlist.ReorderNodes();
	ReorderedNet.m_net = ReorderedNet.ConstructXMG();
	ReorderedNet.ConfigMF();
	ofCSV << ReorderedNet.m_nSize << ", " << ReorderedNet.m_nMF << "," << ReorderedNet.m_nCross << "\n";

	vector<NetList> vecParetoList;
	vecParetoList.push_back(ReorderedNet);

	clock_t beginn = clock();
	cout << "Begin resub\n";
	int nNumSubb = 0;
	NetList netlistt;
	NetList EmptyListt;
	netlistt.m_net = ReorderedNet.m_net.clone();
	netlistt.ConfigWithXMG();
	netlistt.m_vecnSchedule.clear();
	for (int i = 0; i < netlistt.m_vecNode.size(); i++)
		netlistt.m_vecnSchedule.push_back(i);
	netlistt.ConfigMF();
	//cout << netlistt.m_nSize << " " << netlistt.m_nMF << " --\n";
	xmg_network orig_nett = netlistt.m_net.clone();
	while (true)
	{
		if (!ResubMult(netlistt))
			break;
		nNumSubb++;
		//cout << "===== root = " << m_nRoot << "\n";
		//cout << "===== valid fail = " << m_nVFail << "\n";
		cout << nNumSubb << "-th resub\n";
		xmg_network new_nett = netlistt.ConstructXMG();
		netlistt = EmptyListt;
		netlistt.m_net = new_nett.clone();
		netlistt.ConfigWithXMG();
		netlistt.m_vecnSchedule.clear();
		for (int i = 0; i < netlistt.m_vecNode.size(); i++)
			netlistt.m_vecnSchedule.push_back(i);
		netlistt.ConfigMF();
		cout << "New netlist: #nodes = " << netlistt.m_nSize << ", MF = " << netlistt.m_nMF << "\n";
		ofCSV << netlistt.m_nSize << ", " << netlistt.m_nMF << "," << netlistt.m_nCross << "," << 1 << "\n";
		if (UpdatePareto(netlistt, vecParetoList))
			cout << "Into Pareto\n";
		else
			cout << "Not Pareto\n";

		if (netlistt.DelRedun())
		{
			cout << "delete reduntant---------------------\n";
			xmg_network newnew_nett = netlistt.ConstructXMG();
			netlistt = EmptyListt;
			netlistt.m_net = newnew_nett.clone();
			netlistt.ConfigWithXMG();
			netlistt.m_vecnSchedule.clear();
			for (int i = 0; i < netlistt.m_vecNode.size(); i++)
				netlistt.m_vecnSchedule.push_back(i);
			netlistt.ConfigMF();
			cout << "New netlist: #nodes = " << netlistt.m_nSize << ", MF = " << netlistt.m_nMF << "\n";
			ofCSV << netlistt.m_nSize << ", " << netlistt.m_nMF << "," << netlistt.m_nCross << "," << 1 << "\n";
			if (UpdatePareto(netlistt, vecParetoList))
				cout << "Into Pareto\n";
			else
				cout << "Not Pareto\n";
		}
		else
			cout << "no reduntant nodes--------------------\n";
	}
	cout << "End resub\n";
	cout << "Time = " << double(clock() - beginn) / CLOCKS_PER_SEC << ", #resub = " << nNumSubb << "\n";


	for (int i = 1; i < m_nRun; i++)
	{
		cout << "-----------------Round " << i << ":-----------------\n";
		int nPick = rand() % vecParetoList.size();
		NetList& OrigNet = vecParetoList[nPick];
		cout << "Pick " << nPick << "-th pareto design\n";
		int nBegin = OrigNet.m_vecMaxMFIndex[0];
		int nEnd = OrigNet.m_vecMaxMFIndex.back();
		int nCritical = (double)OrigNet.m_nMF * m_critical;
		//cout << "Begin/end " << nBegin << "/" << nEnd << "\n";
		//cout << "critical " << nCritical << "\n";
		while (nBegin > 0)
		{
			nBegin--;
			if (OrigNet.m_vecMF[nBegin] < nCritical)
				break;
		}
		while (nEnd < OrigNet.m_vecNode.size() - 1)
		{
			nEnd++;
			if (OrigNet.m_vecMF[nEnd] < nCritical)
				break;
		}

		//cout << "Begin/end " << nBegin << "/" << nEnd << "\n";
		int nOtherNode = OrigNet.m_vecNode.size() - (nEnd - nBegin + 1);

		vector<int> vecTIOrigIndex;
		vector<int> vecPOOrigIndex;
		NetList Sub = OrigNet.ExtractSub(nBegin, nEnd, vecTIOrigIndex, vecPOOrigIndex);
		xmg_network NewNet = Sub.m_net.clone();
		cout << "Extract " << nEnd - nBegin + 1 << " nodes " << nBegin << " - " << nEnd << " / " << OrigNet.m_vecNode.size() << "\n";
		while (true)
		{
			//int nNodeMin = INT_MAX;
			//xmg_network KeepNet = Sub.m_net.clone();
			for (int j = 0; j < 10; j++)
			{
				int nCommand = rand() % 10;
				//cout << nCommand << " ";
				if (nCommand == 0)
					xmg_depth_rewrite(NewNet, true, 'a', 1.2);
				else if (nCommand == 1)
					xmg_depth_rewrite(NewNet, false, 's');
				else if (nCommand == 2)
					xmg_resub(NewNet);
				else if (nCommand == 3)
					xmg_node_resynthesis(NewNet, 3);
				else if (nCommand == 4)
					xmg_node_resynthesis(NewNet, 4);
				else if (nCommand == 5)
					xmg_cut_rewrite(NewNet, 2);
				else if (nCommand == 6)
					xmg_cut_rewrite(NewNet, 3);
				else if (nCommand == 7)
					xmg_cut_rewrite(NewNet, 4);
				else if (nCommand == 8)
					xmg_resub(NewNet, 5);
				else
				{
					functional_reduction(NewNet);
					NewNet = cleanup_dangling(NewNet);
				}
				/*
				if (NewNet.num_gates() < nNodeMin)
				{
					KeepNet = NewNet.clone();
					nNodeMin = NewNet.num_gates();
				}
				*/
				//if (NewNet.num_gates() > KeepNet.num_gates())
					//NewNet = KeepNet.clone();
				//else
					//KeepNet = NewNet.clone();
				cout << NewNet.num_gates() << "; ";
			}
			int nCh = NewNet.num_gates();
			while (true)
			{
				Resyn2(NewNet);
				if (NewNet.num_gates() >= nCh)
					break;
				nCh = NewNet.num_gates();
			}
			cout << NewNet.num_gates() << "; ";
			cout << "\n";
			if ((NewNet.num_gates() + nOtherNode) <= nAcceptNode)
				break;
			NewNet = Sub.m_net.clone();
		}
		//write_verilog(NewNet, "err.v");
		Scheduler SubScheduler;
		SubScheduler.m_netlist.m_net = NewNet.clone();
		SubScheduler.m_netlist.ConfigWithXMG(nOrigOffset);
		int nNewNodes = NewNet.num_gates() + nOtherNode;
		int nBound = m_nBound;
		for (NetList& nl : vecParetoList)
		{
			if (nl.m_nSize < nNewNodes)
				nBound = min(nBound, nl.m_nMF - 1);
			else if (nl.m_nSize == nNewNodes)
				nBound = min(nBound, nl.m_nMF);
		}
		//cout << "target MF: " << nBound << "\n";
		SubScheduler.m_nBound = nBound * 1.05;//tmp
		//SubScheduler.m_nBound = m_nBound;
		SubScheduler.ThreadIterPartScheduler();
		if (SubScheduler.m_nFootPrint < 1 || SubScheduler.m_bStop)
		{
			cout << "Not a pareto, return\n";
			continue;
		}
		//cout << "Schedule done\n";
		NetList NewSub = SubScheduler.m_netlist.ReorderNodes();
		NetList NewNetList = OrigNet.SubstituteSub(nBegin, nEnd, NewSub, vecTIOrigIndex, vecPOOrigIndex);
		//cout << "Substitute done\n";
		NewNetList.ConfigMF();
		int nNewMF = NewNetList.m_nMF;
		int nNewSize = NewNetList.m_nSize;
		cout << "New netlist: #nodes = " << nNewSize << ", MF = " << nNewMF << "\n";
		ofCSV << nNewSize << ", " << nNewMF << "," << NewNetList.m_nCross << "\n";
		if (UpdatePareto(NewNetList, vecParetoList))
			cout << "Into Pareto\n";
		else
			cout << "Not Pareto\n";

		NetList EmptyList;
		clock_t begin = clock();
		cout << "Begin resub\n";
		int nNumSub = 0;
		NetList netlist = EmptyList;
		netlist.m_net = NewNetList.m_net.clone();
		netlist.ConfigWithXMG();
		netlist.m_vecnSchedule.clear();
		for (int i = 0; i < netlist.m_vecNode.size(); i++)
			netlist.m_vecnSchedule.push_back(i);
		netlist.ConfigMF();
		xmg_network orig_net = netlist.m_net.clone();
		while (true)
		{
			if (!ResubMult(netlist))
				break;
			nNumSubb++;
			nNumSub++;
			cout << nNumSub << "-th resub\n";
			xmg_network new_net = netlist.ConstructXMG();
			netlist = EmptyList;
			netlist.m_net = new_net.clone();
			netlist.ConfigWithXMG();
			netlist.m_vecnSchedule.clear();
			for (int i = 0; i < netlist.m_vecNode.size(); i++)
				netlist.m_vecnSchedule.push_back(i);
			netlist.ConfigMF();
			cout << "New netlist: #nodes = " << netlist.m_nSize << ", MF = " << netlist.m_nMF << "\n";
			ofCSV << netlist.m_nSize << ", " << netlist.m_nMF << "," << netlist.m_nCross << "," << 1 << "\n";
			if (UpdatePareto(netlist, vecParetoList))
				cout << "Into Pareto\n";
			else
				cout << "Not Pareto\n";

			if (netlist.DelRedun())
			{
				cout << "delete reduntant---------------------\n";
				xmg_network newnew_net = netlist.ConstructXMG();
				netlist = EmptyList;
				netlist.m_net = newnew_net.clone();
				netlist.ConfigWithXMG();
				netlist.m_vecnSchedule.clear();
				for (int i = 0; i < netlist.m_vecNode.size(); i++)
					netlist.m_vecnSchedule.push_back(i);
				netlist.ConfigMF();
				cout << "New netlist: #nodes = " << netlist.m_nSize << ", MF = " << netlist.m_nMF << "\n";
				ofCSV << netlist.m_nSize << ", " << netlist.m_nMF << "," << netlist.m_nCross << "," << 1 << "\n";
				if (UpdatePareto(netlist, vecParetoList))
					cout << "Into Pareto\n";
				else
					cout << "Not Pareto\n";
			}
			else
				cout << "no reduntant nodes--------------------\n";
		}
		cout << "End resub\n";
		cout << "Time = " << double(clock() - begin) / CLOCKS_PER_SEC << ", #resub = " << nNumSub << "\n";
		ofCSV.flush();
	}
	ofCSV.close();
}