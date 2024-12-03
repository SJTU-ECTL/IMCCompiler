#include "Scheduler.h"
int main()
{
	ifstream fin;
	fin.open("benchmarks.txt", ios::in);
	string strBench, strLine;
	int nBound;
	istringstream sstream;
	while (getline(fin, strLine))
	{
		nBound = -1;
		sstream.clear();
		sstream.str(strLine);
		sstream >> strBench >> nBound;
		Scheduler MyScheduler;
		MyScheduler.m_netlist.m_strBench = strBench;
		MyScheduler.m_netlist.ReadFromFile(strBench);  //XMG netlist (.v / .bliff / .aig)
		if (nBound == -1)
			MyScheduler.m_nBound = MyScheduler.m_netlist.m_vecNode.size();
		else
			MyScheduler.m_nBound = nBound;
		MyScheduler.NewDSE();
	}
	getchar();
	fin.close();
	return 0;
}