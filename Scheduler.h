#pragma once
#include "NetList.h"

class Scheduler
{
public:
	Scheduler();
	int CallSMT(NetList& netlist, int nMF);
	int BiSMT(NetList& netlist, int nUpper, int nLower);
	int PartitionILP(NetList& netlist, double dPOWeight = 1);
	void IterPart(NetList& netlist, double dLeft = 1, double dRight = 1);
	void ScheduleThread();
	void ThreadIterPartScheduler();

	bool IsDivisor(const Node& ndRoot, const Node& ndDiv);
	void ResubThread();
	bool ResubMult(NetList& netlist);
	bool UpdatePareto(NetList& NewNetList, vector<NetList>& vecParetoList);
	void NewDSE();

	NetList m_netlist;
	int m_nFootPrint;
	vector<NetList> m_vecPartNet;
	vector<int> m_vecMF;
	int m_nMFLow;
	bool m_bStop;

	bool m_bFoundResub;
	volatile LONG m_nCheckPeak;
	vector<kitty::partial_truth_table> m_tts;
	vector<bool> m_vecPhase;
	vector<int> m_vecPeakIndex;
	volatile LONG m_nVFail;
	volatile LONG m_nRoot;
	NetList m_curNetlist;
	NetList m_newNetlist;
	
	CRITICAL_SECTION m_rCritical;
	HANDLE	m_hEventKillThread;
	int m_nTotNet;
	volatile LONG m_nProcNet;

	int m_nBound;
	double m_epsilon;
	int m_nThread;
	int m_nGraphBound;
	int m_nProgressDone;
	int m_nProgressTotal;
	double m_dPO;

	double m_increase;
	double m_critical;
	int m_nRun;
};