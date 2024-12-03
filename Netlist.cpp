#include "Netlist.h"

NetList::NetList(vector<Node> vecNode, vector<Node> vecIn, vector<int> vecnPo)
{
	m_nArrayRow = 253;
	m_nNumArray = 4;
	m_vecNode = vecNode;
	m_vecIn = vecIn;
	m_vecnPO = vecnPo;
}

void NetList::ReadFromFile(string strFile)
{//config XMG from input strFile
	if (strFile != "")
	{
		if (strFile.back() == 'v')
			lorina::read_verilog(strFile, verilog_reader(m_net));
		else if (strFile.back() == 'f')
			m_net = Bliff2Xmg(strFile);
		else
			m_net = Aig2Xmg(strFile);
	}
	ConfigWithXMG();
}

void NetList::ConfigWithXMG(int nOrigOffset)
{//config netlist with XMG
	m_vecIn.clear();
	m_vecNode.clear();
	m_vecnPO.clear();
	m_vecnSchedule.clear();
	unsigned int nNumGates = m_net.num_gates();
	m_nOffset = m_net.size() - nNumGates;
	m_nNumPI = m_nOffset - 1;

	int nNumTI = 0;
	if (m_nOffset > nOrigOffset)
	{
		nNumTI = m_nOffset - nOrigOffset;
		vector<Node> vecNode(nNumTI);
		m_vecIn = vecNode;
		for (unsigned int i = 0; i < nNumTI; i++)
		{
			m_vecIn[i].m_nIndex = i + nNumGates;
			m_vecIn[i].m_nOrigIndex = i + nNumGates;
		}
	}

	vector<Node> vecNode(nNumGates);
	m_vecNode = vecNode;
	for (unsigned int i = 0; i < nNumGates; i++)
	{
		m_vecNode[i].m_nIndex = i;
		m_vecNode[i].m_nOrigIndex = i;
	}

	m_net.foreach_gate([&](auto const& n)
		{
			int nIndexNow = m_net.node_to_index(n) - m_nOffset;
			if (m_net.is_maj(n))
				m_vecNode[nIndexNow].m_bMAJ = true;
			else if (m_net.is_xor3(n))
				m_vecNode[nIndexNow].m_bMAJ = false;
			else
				cout << "ERROR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";

			m_net.foreach_fanin(n, [&](auto const& f)
				{
					bool bCompf = m_net.is_complemented(f);
					auto fNode = m_net.get_node(f);
					int fIndex = m_net.node_to_index(fNode);
					//constant or PI or TI or regular
					if (m_net.is_constant(fNode))//constant
					{
						if (bCompf)
							m_vecNode[nIndexNow].m_nConstPI = 1;
						else
							m_vecNode[nIndexNow].m_nConstPI = 0;
					}//PI or TI or regular
					else if (!m_net.is_pi(fNode))//regular
					{
						int fIndexNow = fIndex - m_nOffset;
						m_vecNode[fIndexNow].m_vecnSucc.push_back(nIndexNow);
						m_vecNode[nIndexNow].m_vecnPred.push_back(fIndexNow);
						//m_vecNode[fIndexNow].m_vecNodeSucc.push_back(&m_vecNode[nIndexNow]);
						//m_vecNode[nIndexNow].m_vecNodePred.push_back(&m_vecNode[fIndexNow]);
						m_vecNode[nIndexNow].m_vecbPredComp.push_back(bCompf);
					}//PI or TI
					else if (fIndex >= nOrigOffset) //TI
					{
						m_vecIn[fIndex - nOrigOffset].m_vecnSucc.push_back(nIndexNow);
						m_vecNode[nIndexNow].m_vecnPred.push_back(fIndex - nOrigOffset + nNumGates);
						m_vecNode[nIndexNow].m_vecbPredComp.push_back(bCompf);
					}//PI
					else
					{
						m_vecNode[nIndexNow].m_vecnPredPI.push_back(fIndex - 1);
						m_vecNode[nIndexNow].m_vecbPredPIComp.push_back(bCompf);
					}
				});
		});

	m_net.foreach_po([&](auto const& f, auto i)
		{
			int nPO = m_net.node_to_index(m_net.get_node(f)) - m_nOffset;
			int nTI = m_net.node_to_index(m_net.get_node(f)) - nOrigOffset;
			if ((nPO >= 0) && !IsInVector(nPO, m_vecnPO))
			{
				m_vecnPO.push_back(nPO);
				m_vecNode[nPO].m_bPO = true;
			}
			if ((nPO < 0) && (nTI >= 0) && !IsInVector(nTI + nNumGates, m_vecnPO))
			{
				m_vecnPO.push_back(nTI + nNumGates);
				m_vecIn[nTI].m_bPO = true;
			}
			m_vecnPOIndex.push_back(m_net.node_to_index(m_net.get_node(f)) - 1);
			m_vecbPOComp.push_back(m_net.is_complemented(f));
		});
	if (nOrigOffset == 1000000)
	{
		for (int i = 0; i < m_vecNode.size(); i++)
		{
			Node& nd = m_vecNode[i];
			for (int pi : nd.m_vecnPredPI)
				nd.m_setnConePI.insert(pi);
			for (int fi : nd.m_vecnPred)
				nd.m_setnConePI.insert(m_vecNode[fi].m_setnConePI.begin(), m_vecNode[fi].m_setnConePI.end());
		}
	}
}

void NetList::ConfigMF()
{//config MF with scheduling result
	m_vecMF.clear();
	int nNumNode = m_vecNode.size();
	m_nSize = nNumNode;
	m_nNumPI = m_nOffset - 1;
	vector<int> vecNumUnSchedSuc(nNumNode, 0);
	vector<int> vecMemStatus;
	for (unsigned int n = 0; n < nNumNode; n++)
		vecNumUnSchedSuc[n] = m_vecNode[n].m_vecnSucc.size();
	m_vecMaxMFIndex.clear();
	int nMaxMF = -1;
	int nMemRow;
	int nCurrentMF = 0;
	m_nCross = 0;
	for (int& n : m_vecnSchedule)
	{
		Node& node = m_vecNode[n];
		vector<int> vecPredArray;
		for (int& i : node.m_vecnPred)
		{
			int nPredMem = GetIndexInVector(i, vecMemStatus);
			if (nPredMem == -1)
			{
				m_nMF = -1;
				return;
			}
			vecPredArray.push_back((int)(nPredMem + m_nOffset) / m_nArrayRow);
			vecNumUnSchedSuc[i]--;
			if (vecNumUnSchedSuc[i] == 0 && !m_vecNode[i].m_bPO)
			{
				vecMemStatus[nPredMem] = -1;
				nCurrentMF--;
			}
		}

		for (nMemRow = 0; nMemRow < vecMemStatus.size(); nMemRow++)
		{
			if (vecMemStatus[nMemRow] == -1)
				break;
		}
		
		if (nMemRow == vecMemStatus.size())
			vecMemStatus.push_back(-1);
		vecMemStatus[nMemRow] = n;
		int nNowArray = (int)(nMemRow + m_nOffset) / m_nArrayRow;

		for (int i : vecPredArray)
		{
			if (i != nNowArray)
				m_nCross++;
		}

		for (int nIn : node.m_vecnPredPI)
		{
			if (nNowArray != ((int)nIn / m_nArrayRow))
				m_nCross++;
		}

		nCurrentMF++;
		m_vecMF.push_back(nCurrentMF);
		if (nCurrentMF < nMaxMF)
			continue;
		if (nCurrentMF > nMaxMF)
		{
			m_vecMaxMFIndex.clear();
			nMaxMF = nCurrentMF;
		}
		m_vecMaxMFIndex.push_back(m_vecMF.size() - 1);

	}
	for (int i : m_vecnPO)
	{
		int nMem = GetIndexInVector(i, vecMemStatus);
		if (nMem == -1)
		{
			m_nMF = -1;
			return;
		}
	}
	m_nMF = vecMemStatus.size();
	cout << "Checking result: Size = " << m_nSize << "; MF = " << m_nMF << "\n";
}

xmg_network NetList::ConstructXMG(int nOrigOffset)
{//construct XMG with netlist
	int nNumNode = m_vecNode.size();
	int nNumPI = (m_nOffset > nOrigOffset) ? nOrigOffset - 1 : m_nOffset - 1;
	int nNumTI = (m_nOffset > nOrigOffset) ? m_nOffset - nOrigOffset : 0;
	int nNumIn = m_nOffset - 1;
	xmg_network ntk;
	xmg_network::signal sigConst = ntk.get_constant(false);
	vector<xmg_network::signal> vecPI;
	for (int i = 0; i < nNumPI; i++)
		vecPI.push_back(ntk.create_pi());
	vector<xmg_network::signal> vecTI;
	for (int i = 0; i < nNumTI; i++)
		vecTI.push_back(ntk.create_pi());
	vector<xmg_network::signal> vecNode;
	for (Node& nd : m_vecNode)
	{
		//cout << "-node " << nd.m_nIndex << "\n";
		vector<xmg_network::signal> vecFI;
		if (nd.m_nConstPI == 1)
			vecFI.push_back(ntk.create_not(sigConst));
		if (nd.m_nConstPI == 0)
			vecFI.push_back(sigConst);
		for (int i = 0; i < nd.m_vecnPredPI.size(); i++)
		{
			if (nd.m_vecbPredPIComp[i])
				vecFI.push_back(ntk.create_not(vecPI[nd.m_vecnPredPI[i]]));
			else
				vecFI.push_back(vecPI[nd.m_vecnPredPI[i]]);
		}
		for (int i = 0; i < nd.m_vecnPred.size(); i++)
		{
			if (nd.m_vecnPred[i] >= nNumNode)//TI
			{
				int nTIIndex = nd.m_vecnPred[i] - nNumNode;
				if (nd.m_vecbPredComp[i])
					vecFI.push_back(ntk.create_not(vecTI[nTIIndex]));
				else
					vecFI.push_back(vecTI[nTIIndex]);
			}
			else
			{
				if (nd.m_vecbPredComp[i])
					vecFI.push_back(ntk.create_not(vecNode[nd.m_vecnPred[i]]));
				else
					vecFI.push_back(vecNode[nd.m_vecnPred[i]]);
			}
		}
		if (vecFI.size() != 3)
			cout << "ERROR!!!!SIZE!!!!!!!!!!\n";
		if (nd.m_bMAJ)
			vecNode.push_back(ntk.create_maj(vecFI[0], vecFI[1], vecFI[2]));
		else
			vecNode.push_back(ntk.create_xor3(vecFI[0], vecFI[1], vecFI[2]));
		//cout << ntk.num_gates() << "\n";
	}
	for (int i = 0; i < m_vecnPOIndex.size(); i++)
	{
		int nIndex = m_vecnPOIndex[i];
		if (nIndex == -1)
		{
			ntk.create_po(m_vecbPOComp[i] ? ntk.create_not(sigConst) : sigConst);
			continue;
		}
		if (nIndex < nNumPI)
			ntk.create_po(m_vecbPOComp[i] ? ntk.create_not(vecPI[nIndex]) : vecPI[nIndex]);
		else if (nIndex < (nNumPI + nNumTI))
			ntk.create_po(m_vecbPOComp[i] ? ntk.create_not(vecTI[nIndex - nNumPI]) : vecTI[nIndex - nNumPI]);
		else
			ntk.create_po(m_vecbPOComp[i] ? ntk.create_not(vecNode[nIndex - nNumIn]) : vecNode[nIndex - nNumIn]);
	}
	return ntk;
}

NetList NetList::ReorderNodes()
{//reorder nodes according to scheduling result
	int nNumNode = m_vecNode.size();
	int nNumIn = m_vecIn.size();
	if (m_vecnSchedule.size() != nNumNode)
	{
		cout << "No Schedule\n";
		return NetList();
	}
	vector<Node> vecNode;
	vector<Node> vecIn;
	vector<int> vecnPo = m_vecnPO;
	std::map<int, int> mapOldNew;
	std::map<int, int> mapOldNewPOIndex;
	for (int i = 0; i < nNumNode; i++)
	{
		mapOldNew[m_vecnSchedule[i]] = i;
		mapOldNewPOIndex[m_vecnSchedule[i] + m_nOffset - 1] = i + m_nOffset - 1;
	}
	for (Node ndOld : m_vecIn)
	{
		Node ndNew = ndOld;
		RemapVector(ndNew.m_vecnSucc, mapOldNew);
		vecIn.push_back(ndNew);
	}
	for (int i = 0; i < nNumNode; i++)
	{
		Node ndNew = m_vecNode[m_vecnSchedule[i]];
		RemapVector(ndNew.m_vecnSucc, mapOldNew);
		RemapVector(ndNew.m_vecnPred, mapOldNew);
		ndNew.m_nIndex = i;
		ndNew.m_nOrigIndex = i;
		vecNode.push_back(ndNew);
	}
	RemapVector(vecnPo, mapOldNew);
	NetList NewNetList(vecNode, vecIn, vecnPo);
	NewNetList.m_nOffset = m_nOffset;
	for (int i = 0; i < nNumNode; i++)
		NewNetList.m_vecnSchedule.push_back(i);
	NewNetList.m_vecnPOIndex = m_vecnPOIndex;
	RemapVector(NewNetList.m_vecnPOIndex, mapOldNewPOIndex);
	NewNetList.m_vecbPOComp = m_vecbPOComp;
	//NewNetList.Print();
	//NewNetList.m_net = NewNetList.ConstructXMG();
	//cout << "reorder done\n";
	//cout << NewNetList.m_vecNode.size() << "\n";
	return NewNetList;
}

vector<int> NetList::GetPeak()
{//get peak in netlist
	vector<int> vecPeak;
	int nNumNode = m_vecNode.size();
	vector<int> vecNumUnSchedSuc(nNumNode, 0);
	vector<int> vecMemStatus;
	for (unsigned int n = 0; n < nNumNode; n++)
		vecNumUnSchedSuc[n] = m_vecNode[n].m_vecnSucc.size();
	int nMaxMF = -1;
	int nMemRow;
	int nCurrentMF = 0;
	for (int& n : m_vecnSchedule)
	{
		Node& node = m_vecNode[n];
		for (int& i : node.m_vecnPred)
		{
			int nPredMem = GetIndexInVector(i, vecMemStatus);
			vecNumUnSchedSuc[i]--;
			if (vecNumUnSchedSuc[i] == 0 && !m_vecNode[i].m_bPO)
			{
				vecMemStatus[nPredMem] = -1;
				nCurrentMF--;
			}
		}
		for (nMemRow = 0; nMemRow < vecMemStatus.size(); nMemRow++)
		{
			if (vecMemStatus[nMemRow] == -1)
				break;
		}
		if (nMemRow == vecMemStatus.size())
			vecMemStatus.push_back(-1);
		vecMemStatus[nMemRow] = n;
		nCurrentMF++;
		if (nCurrentMF > nMaxMF)
		{
			vecPeak = vecMemStatus;
			nMaxMF = nCurrentMF;
		}
	}
	return vecPeak;
}

void NetList::SubPred(int nNode, vector<int>& vecnPredID, vector<bool>& vecbPredComp, bool bMaj, bool bComp)
{//substitute the predecessor for nNode
	//cout << "Node" << nNode << "sub with ID: " << vecnPredID[0] << " " << vecnPredID[1] << " " << vecnPredID[2] << "\n";
	Node& nd = m_vecNode[nNode];
	for (int pred : nd.m_vecnPred)
		m_vecNode[pred].m_vecnSucc.erase(find(m_vecNode[pred].m_vecnSucc.begin(), m_vecNode[pred].m_vecnSucc.end(), nNode));

	nd.m_bMAJ = bMaj;
	nd.m_nConstPI = -1;
	nd.m_vecnPredPI.clear();
	nd.m_vecbPredPIComp.clear();
	nd.m_vecnPred.clear();
	nd.m_vecbPredComp.clear();


	for (int i = 0; i < 3; i++)
	{
		if (vecnPredID[i] >= m_nOffset)
		{
			nd.m_vecnPred.push_back(vecnPredID[i] - m_nOffset);
			nd.m_vecbPredComp.push_back(vecbPredComp[i]);
		}
		else if (vecnPredID[i] > 0)
		{
			nd.m_vecnPredPI.push_back(vecnPredID[i] - 1);
			nd.m_vecbPredPIComp.push_back(vecbPredComp[i]);
		}
		else if (vecbPredComp[i])
			nd.m_nConstPI = 1;
		else
			nd.m_nConstPI = 0;
	}
	for (int pred : nd.m_vecnPred)
		m_vecNode[pred].m_vecnSucc.push_back(nNode);
	if (bComp)
	{
		for (int suc : nd.m_vecnSucc)
		{
			Node& nsuc = m_vecNode[suc];
			for (int i = 0; i < nsuc.m_vecnPred.size(); i++)
			{
				if (nsuc.m_vecnPred[i] == nNode)
					nsuc.m_vecbPredComp[i] = !nsuc.m_vecbPredComp[i];
			}
		}
		if (nd.m_bPO)
		{
			for (int i = 0; i < m_vecnPOIndex.size(); i++)
			{
				if (m_vecnPOIndex[i] == nNode + m_nOffset - 1)
					m_vecbPOComp[i] = !m_vecbPOComp[i];
			}
		}
	}
}

void NetList::SubNode(int nSmall, int nLarge, bool bComp)
{//substitute nLarge with nSmall
	//cout << "SubNode " << nSmall << " " << nLarge << " " << bComp << "\n";
	for (Node& nd : m_vecNode)
	{
		int ndel = -1;
		for (int i = 0; i < nd.m_vecnSucc.size(); i++)
		{
			if (nd.m_vecnSucc[i] == nLarge)
				ndel = i;
			else if (nd.m_vecnSucc[i] > nLarge)
				nd.m_vecnSucc[i]--;
		}
		if (ndel != -1)
			nd.m_vecnSucc.erase(nd.m_vecnSucc.begin() + ndel);
		for (int i = 0; i < nd.m_vecnPred.size(); i++)
		{
			if (nd.m_vecnPred[i] == nLarge)
			{
				nd.m_vecnPred[i] = nSmall;
				if (bComp)
					nd.m_vecbPredComp[i] = !nd.m_vecbPredComp[i];
			}
			else if (nd.m_vecnPred[i] > nLarge)
				nd.m_vecnPred[i]--;
		}
		if (nd.m_nIndex > nLarge)
			nd.m_nIndex--;
		nd.m_nOrigIndex = nd.m_nIndex;
	}

	if (m_vecNode[nLarge].m_bPO && m_vecNode[nSmall].m_bPO)
		m_vecnPO.erase(find(m_vecnPO.begin(), m_vecnPO.end(), nLarge));
	for (int i = 0; i < m_vecnPO.size(); i++)
	{
		if (m_vecnPO[i] == nLarge)
			m_vecnPO[i] = nSmall;
		else if (m_vecnPO[i] > nLarge)
			m_vecnPO[i]--;
	}
	for (int i = 0; i < m_vecnPOIndex.size(); i++)
	{
		if (m_vecnPOIndex[i] == nLarge + m_nOffset - 1)
		{
			m_vecnPOIndex[i] = nSmall + m_nOffset - 1;
			if (bComp)
				m_vecbPOComp[i] = !m_vecbPOComp[i];
		}
		else if (m_vecnPOIndex[i] > nLarge + m_nOffset - 1)
			m_vecnPOIndex[i]--;
	}

	m_vecNode.erase(m_vecNode.begin() + nLarge);
}

bool NetList::DelRedun()
{//delete redunant
	bool bChange = false;
	while (true)
	{
		int nRed = 0;
		for (nRed = 0; nRed < m_vecNode.size(); nRed++)
		{
			if ((!m_vecNode[nRed].m_bPO) && (m_vecNode[nRed].m_vecnSucc.size() == 0))
				break;
		}
		if (nRed == m_vecNode.size())
			break;
		//cout << "delete node " << nRed << "\n";
		bChange = true;
		for (Node& nd : m_vecNode)
		{
			int ndel = -1;
			for (int i = 0; i < nd.m_vecnSucc.size(); i++)
			{
				if (nd.m_vecnSucc[i] == nRed)
					ndel = i;
				else if (nd.m_vecnSucc[i] > nRed)
					nd.m_vecnSucc[i]--;
			}
			if (ndel != -1)
				nd.m_vecnSucc.erase(nd.m_vecnSucc.begin() + ndel);
			for (int i = 0; i < nd.m_vecnPred.size(); i++)
			{
				if (nd.m_vecnPred[i] == nRed)
					cout << "fanin still exist!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
				else if (nd.m_vecnPred[i] > nRed)
					nd.m_vecnPred[i]--;
			}
			if (nd.m_nIndex > nRed)
				nd.m_nIndex--;
			nd.m_nOrigIndex = nd.m_nIndex;
		}
		m_vecNode.erase(m_vecNode.begin() + nRed);
		for (int i = 0; i < m_vecnPO.size(); i++)
		{
			if (m_vecnPO[i] > nRed)
				m_vecnPO[i]--;
		}
		for (int i = 0; i < m_vecnPOIndex.size(); i++)
		{
			if (m_vecnPOIndex[i] > nRed + m_nOffset - 1)
				m_vecnPOIndex[i]--;
		}
	}
	return bChange;
}


NetList NetList::ExtractSub(int nBegin, int nEnd, vector<int>& vecTIOrigIndex, vector<int>& vecPOOrigIndex)
{//extract sub-netlist
	vecTIOrigIndex.clear();
	vecPOOrigIndex.clear();
	vector<Node> vecNode;
	vector<Node> vecIn;
	vector<int> vecnPO;
	vector<int> vecnPOIndex;
	int nOrigNodeNum = m_vecNode.size();
	int nSubNodeNum = nEnd - nBegin + 1;
	vector<int> vecNumUnSchedSuc(nEnd + 1, 0);
	for (unsigned int n = 0; n <= nEnd; n++)
		vecNumUnSchedSuc[n] = m_vecNode[n].m_vecnSucc.size();
	set<int> setMemory;
	//build map, original index to index in subgraph
	std::map<int, int> mapOrigSub;
	for (int n = 0; n < nSubNodeNum; n++)
		mapOrigSub[n + nBegin] = n;
	for (int n = 0; n <= nEnd; n++)
	{
		if (n == nBegin)
		{
			int nTIIndex = nSubNodeNum;//TI index begin from #nodes
			for (const int& i : setMemory)
			{
				mapOrigSub[i] = nTIIndex;
				vecTIOrigIndex.push_back(i);
				nTIIndex++;
			}
		}
		Node& node = m_vecNode[n];
		for (int& i : node.m_vecnPred)
		{
			vecNumUnSchedSuc[i]--;
			if (vecNumUnSchedSuc[i] == 0 && !m_vecNode[i].m_bPO)
				setMemory.erase(i);
		}
		setMemory.insert(n);
	}
	for (const int& nOrigIndex : vecTIOrigIndex)
	{
		int nSubIndex = mapOrigSub[nOrigIndex];
		Node ndNew = m_vecNode[nOrigIndex];
		RemapVector(ndNew.m_vecnSucc, mapOrigSub, true);
		ndNew.m_vecnPred.clear();
		ndNew.m_vecbPredComp.clear();
		ndNew.m_vecnPredPI.clear();
		ndNew.m_vecbPredPIComp.clear();
		ndNew.m_nIndex = nSubIndex;
		ndNew.m_nOrigIndex = nSubIndex;
		if (setMemory.find(nOrigIndex) != setMemory.end())
		{
			ndNew.m_bPO = true;
			vecnPO.push_back(nSubIndex);
			vecnPOIndex.push_back(nSubIndex - nSubNodeNum + m_nOffset - 1);
			vecPOOrigIndex.push_back(nOrigIndex);
		}
		else
			ndNew.m_bPO = false;
		vecIn.push_back(ndNew);
	}
	for (int nOrigIndex = nBegin; nOrigIndex <= nEnd; nOrigIndex++)
	{
		int nSubIndex = mapOrigSub[nOrigIndex];
		Node ndNew = m_vecNode[nOrigIndex];
		RemapVector(ndNew.m_vecnSucc, mapOrigSub, true);
		RemapVector(ndNew.m_vecnPred, mapOrigSub);
		ndNew.m_nIndex = nSubIndex;
		ndNew.m_nOrigIndex = nSubIndex;
		if (setMemory.find(nOrigIndex) != setMemory.end())
		{
			ndNew.m_bPO = true;
			vecnPO.push_back(nSubIndex);
			vecnPOIndex.push_back(nSubIndex + m_nOffset - 1 + vecTIOrigIndex.size());
			vecPOOrigIndex.push_back(nOrigIndex);
		}
		else
			ndNew.m_bPO = false;
		vecNode.push_back(ndNew);
	}
	NetList NewNetList(vecNode, vecIn, vecnPO);
	NewNetList.m_nOffset = m_nOffset + vecIn.size();
	NewNetList.m_vecnPOIndex = vecnPOIndex;
	vector<bool> tmp(vecnPOIndex.size(), false);
	NewNetList.m_vecbPOComp = tmp;
	NewNetList.m_net = NewNetList.ConstructXMG(m_nOffset);
	return NewNetList;
}

NetList NetList::SubstituteSub(int nBegin, int nEnd, const NetList& NewSub, vector<int>& vecTIOrigIndex, vector<int>& vecPOOrigIndex)
{//substitute sub-netlist
	//cout << "Substituting from " << nBegin << " to " << nEnd << "with new sub of " << NewSub.m_vecNode.size() << "nodes\n";
	vector<Node> vecNode;
	vector<Node> vecIn;
	vector<int> vecnPO;
	int nNewSubNode = NewSub.m_vecNode.size();
	int nBackOffset = nNewSubNode - (nEnd - nBegin + 1);
	for (int n = 0; n < nBegin; n++)
	{
		Node NewNode = m_vecNode[n];
		vector<int>::iterator it = find(vecTIOrigIndex.begin(), vecTIOrigIndex.end(), n);
		if (it != vecTIOrigIndex.end())//TI in sub graph
		{
			int nTIIndex = distance(vecTIOrigIndex.begin(), it);
			NewNode.m_vecnSucc.clear();
			for (int i : m_vecNode[n].m_vecnSucc)
			{
				if (i < nBegin)
					NewNode.m_vecnSucc.push_back(i);
				else if (i > nEnd)
					NewNode.m_vecnSucc.push_back(i + nBackOffset);
			}
			for (int i : NewSub.m_vecIn[nTIIndex].m_vecnSucc)
			{
				if (i < nNewSubNode)
					NewNode.m_vecnSucc.push_back(i + nBegin);
			}
		}
		vecNode.push_back(NewNode);
	}
	std::map<int, int> mapSubPOIndexOrigNew;
	std::map<int, bool> mapSubPOIndexOrigComp;
	for (int n = 0; n < nNewSubNode; n++)
	{
		Node NewNode = NewSub.m_vecNode[n];
		NewNode.m_nIndex = n + nBegin;
		NewNode.m_nOrigIndex = n + nBegin;
		for (int i = 0; i < NewNode.m_vecnPred.size(); i++)
		{
			if (NewNode.m_vecnPred[i] < nNewSubNode)
				NewNode.m_vecnPred[i] += nBegin;
			else
				NewNode.m_vecnPred[i] = vecTIOrigIndex[NewNode.m_vecnPred[i] - nNewSubNode];
		}
		for (int i = 0; i < NewNode.m_vecnSucc.size(); i++)
			NewNode.m_vecnSucc[i] += nBegin;
		if (NewNode.m_bPO)
		{
			int nSubPOPosition = distance(NewSub.m_vecnPO.begin(), find(NewSub.m_vecnPO.begin(), NewSub.m_vecnPO.end(), n));
			int nOrigIndex = vecPOOrigIndex[nSubPOPosition];
			mapSubPOIndexOrigNew[nOrigIndex] = NewNode.m_nIndex;
			mapSubPOIndexOrigComp[nOrigIndex] = NewSub.m_vecbPOComp[nSubPOPosition];
			Node& oldNode = m_vecNode[nOrigIndex];
			NewNode.m_bPO = oldNode.m_bPO;
			for (int i : oldNode.m_vecnSucc)
			{
				if (i > nEnd)
					NewNode.m_vecnSucc.push_back(i + nBackOffset);
			}
		}
		vecNode.push_back(NewNode);
	}
	for (int n = nEnd + 1; n < m_vecNode.size(); n++)
	{
		Node NewNode = m_vecNode[n];
		NewNode.m_nIndex = n + nBackOffset;
		NewNode.m_nOrigIndex = n + nBackOffset;
		for (int i = 0; i < NewNode.m_vecnPred.size(); i++)
		{
			if (NewNode.m_vecnPred[i] > nEnd)
				NewNode.m_vecnPred[i] += nBackOffset;
			else if (NewNode.m_vecnPred[i] >= nBegin)
			{
				if (mapSubPOIndexOrigComp[NewNode.m_vecnPred[i]])
					NewNode.m_vecbPredComp[i] = !NewNode.m_vecbPredComp[i];
				NewNode.m_vecnPred[i] = mapSubPOIndexOrigNew[NewNode.m_vecnPred[i]];
			}
		}
		for (int i = 0; i < NewNode.m_vecnSucc.size(); i++)
			NewNode.m_vecnSucc[i] += nBackOffset;
		vecNode.push_back(NewNode);
	}
	for (Node& nd : vecNode)
	{
		if (nd.m_bPO)
			vecnPO.push_back(nd.m_nIndex);
	}
	vector<int> vecnPOIndex = m_vecnPOIndex;
	vector<bool> vecbPOComp = m_vecbPOComp;
	for (int n = 0; n < vecnPOIndex.size(); n++)
	{
		if (vecnPOIndex[n] > nEnd + m_nOffset - 1)
			vecnPOIndex[n] += nBackOffset;
		else if (vecnPOIndex[n] >= nBegin + m_nOffset - 1)
		{
			vecnPOIndex[n] = mapSubPOIndexOrigNew[vecnPOIndex[n] - (m_nOffset - 1)] + m_nOffset - 1;
			if (mapSubPOIndexOrigComp[vecnPOIndex[n]])
				vecbPOComp[n] = !vecbPOComp[n];
		}
	}
	NetList NewNetList(vecNode, vecIn, vecnPO);
	NewNetList.m_nOffset = m_nOffset;
	NewNetList.m_vecnPOIndex = vecnPOIndex;
	NewNetList.m_vecbPOComp = vecbPOComp;
	NewNetList.m_net = NewNetList.ConstructXMG();
	NewNetList.m_vecnSchedule.clear();
	for (int i = 0; i < vecNode.size(); i++)
		NewNetList.m_vecnSchedule.push_back(i);
	return NewNetList;
}