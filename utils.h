#pragma once
#include <iostream>
#include <time.h>
#include <vector>
#include <map>
#include <queue>
#include <set>
#include <sstream>
#include<string>
#include<fstream>
#include <cmath>
#include<stack>
#include <algorithm>
#include <random>
#include "mockturtle\mockturtle.hpp"
#include "lorina\aiger.hpp"
#include "C:\\gurobi1001\\win64\\include\\gurobi_c++.h"
#include "z3include\\z3++.h"
#pragma comment(lib, "C:\\z3\\libZ3.lib")
#pragma comment(lib, "gurobi100.lib")
#pragma comment(lib, "gurobi_c++md2017.lib")

using namespace mockturtle;
using namespace std;

inline xmg_network Aig2Xmg(string strAigName)
{
	mockturtle::aig_network aig;
	lorina::read_aiger(strAigName, mockturtle::aiger_reader(aig));
	xmg_npn_resynthesis resyn;
	exact_library<xmg_network, xmg_npn_resynthesis> lib(resyn);
	map_params ps;
	ps.skip_delay_round = true;
	ps.required_time = numeric_limits<double>::max();
	xmg_network xmg = mockturtle::map(aig, lib, ps);

	functional_reduction(xmg);
	xmg = cleanup_dangling(xmg);

	if (strAigName == "bench\\sqrt.aig")
		return xmg;

	depth_view xmg_depth{ xmg };
	fanout_view xmg_fanout{ xmg_depth };
	xmg_resubstitution(xmg_fanout);
	xmg_network xmg_res = xmg_fanout;
	xmg_res = cleanup_dangling(xmg);
	if (xmg_res.size() < xmg.size())
		xmg = xmg_res;

	return xmg;
}

inline xmg_network Bliff2Xmg(string strBliffName)
{
	cover_network cover;
	lorina::read_blif(strBliffName, blif_reader(cover));
	aig_network aig;
	convert_cover_to_graph(aig, cover);

	xmg_npn_resynthesis resyn;
	exact_library<xmg_network, xmg_npn_resynthesis> lib(resyn);

	map_params ps;
	ps.skip_delay_round = true;
	ps.required_time = numeric_limits<double>::max();
	xmg_network xmg = mockturtle::map(aig, lib, ps);
	functional_reduction(xmg);
	xmg = cleanup_dangling(xmg);

	depth_view xmg_depth{ xmg };
	fanout_view xmg_fanout{ xmg_depth };
	xmg_resubstitution(xmg_fanout);
	xmg_network xmg_res = xmg_fanout;
	xmg_res = cleanup_dangling(xmg);
	if (xmg_res.size() < xmg.size())
		xmg = xmg_res;

	return xmg;
}

inline void Trim(string& s)
{//remove ' ' and '\\'
	s.erase(std::remove(s.begin(), s.end(), ' '), s.end());
	s.erase(std::remove(s.begin(), s.end(), '\\'), s.end());
}

inline string ExtractStr(string& s, string strLeft, string strRight)
{//extract substring from strLeft to strRight in s
	size_t left = s.find(strLeft);
	if (left == string::npos)
		return "";
	size_t right = s.find(strRight, left + strLeft.length());
	if (right == string::npos)
		return "";
	//cout << left << " " << right << "\n";
	string res = s.substr(left + strLeft.length(), right - left - strLeft.length());
	s.erase(left, right - left);
	return res;
}

inline int GetIndexInVector(int n, vector<int>& vecn)
{//get index of n in vecn
	for (int i = 0; i < vecn.size(); i++)
	{
		if (vecn[i] == n)
			return i;
	}
	return -1;
}

inline int GetFromArray(vector<int>& vecMem, int nArray, int nArraySize = 256)
{//get empty row in nArray
	int nRow = nArray * nArraySize;
	int nBack = (nArray + 1) * nArraySize;
	for (; nRow < nBack; nRow++)
	{
		if (vecMem[nRow] == 0)
			break;
	}
	if (nRow == nBack)
		nRow = -1;
	return nRow;
}

inline bool IsInVector(int n, vector<int>& vecn)
{//whether n is in vecn
	for (int i = 0; i < vecn.size(); i++)
	{
		if (vecn[i] == n)
			return true;
	}
	return false;
}

inline void RemapVector(vector<int>& vecN, std::map<int, int> mapOldNew, bool bSkip = false)
{//remap vecN with mapOldNew
	for (int i = 0; i < vecN.size(); i++)
	{
		if (mapOldNew.find(vecN[i]) == mapOldNew.end())
		{
			if (bSkip)
				vecN[i] = -100;
			continue;
		}
		vecN[i] = mapOldNew[vecN[i]];
	}
	if (bSkip)
		vecN.erase(remove(vecN.begin(), vecN.end(), -100), vecN.end());
}

inline float xmg_depth_rewrite(xmg_network& xmg, bool allow_size_increase = false, char strat = '0', float overhead = 1.0)
{
	depth_view xmg_depth{ xmg };

	xmg_algebraic_depth_rewriting_params ps;
	cut_rewriting_stats st;
	if (strat == 's')
		ps.strategy = ps.selective;
	else if (strat == 'a')
		ps.strategy = ps.aggressive;
	else
		ps.strategy = ps.dfs;
	ps.overhead = overhead;
	ps.allow_area_increase = allow_size_increase;

	xmg_algebraic_depth_rewriting(xmg_depth, ps);

	xmg = xmg_depth;

	xmg = cleanup_dangling(xmg);
	return to_seconds(st.time_total);
}

inline void xmg_resub(xmg_network& xmg, int max_insert = 2)
{
	depth_view xmg_depth{ xmg };
	fanout_view xmg_fanout{ xmg_depth };
	resubstitution_params ps;
	ps.max_inserts = max_insert;
	xmg_resubstitution(xmg_fanout, ps);
	xmg = xmg_fanout;
	xmg = cleanup_dangling(xmg);
}

inline float xmg_node_resynthesis(xmg_network& xmg, int cut_size = 4)
{
	mapping_view<xmg_network, true> mapped_xmg{ xmg };
	lut_mapping_params ps;
	node_resynthesis_stats st;
	ps.cut_enumeration_ps.cut_size = cut_size;
	lut_mapping<mapping_view<xmg_network, true>, true>(mapped_xmg, ps);

	const auto klut = *collapse_mapped_network<klut_network>(mapped_xmg);
	xmg_npn_resynthesis resyn;
	xmg = node_resynthesis<xmg_network>(klut, resyn);
	xmg = cleanup_dangling(xmg);
	return to_seconds(st.time_total);
}

inline float xmg_cut_rewrite(xmg_network& xmg, int cut_size = 4, bool bZero = true)
{
	exact_xmg_resynthesis<xmg_network> resyn;
	cut_rewriting_params ps;
	cut_rewriting_stats st;
	ps.cut_enumeration_ps.cut_size = cut_size;
	ps.allow_zero_gain = bZero;

	xmg = cut_rewriting(xmg, resyn, ps);
	xmg = cleanup_dangling(xmg);
	return to_seconds(st.time_total);
}

inline void Balance(xmg_network& xmg)
{
	xmg_network n = xmg.clone();
	sop_rebalancing<xmg_network> balance_fn;
	n = balancing(n, { balance_fn });
	if (n.num_gates() <= xmg.num_gates())
		xmg = n.clone();
}

inline void Rewrite(xmg_network& xmg, bool bZero = false)
{
	xmg_network n = xmg.clone();
	exact_xmg_resynthesis<xmg_network> resyn;
	cut_rewriting_params ps;
	ps.cut_enumeration_ps.cut_size = 4;
	ps.allow_zero_gain = bZero;
	n = cut_rewriting(n, resyn, ps);
	n = cleanup_dangling(n);
	if (n.num_gates() <= xmg.num_gates())
		xmg = n.clone();
}

inline void Refactor(xmg_network& xmg, bool bZero = false)
{
	xmg_network n = xmg.clone();
	akers_resynthesis<xmg_network> resyn;
	refactoring_params ps;
	ps.allow_zero_gain = bZero;
	refactoring(n, resyn, ps);
	n = cleanup_dangling(n);
	if (n.num_gates() <= xmg.num_gates())
		xmg = n.clone();
}

inline void Resyn2(xmg_network& xmg)
{
	Balance(xmg);
	//cout << xmg.num_gates() << " ";
	Rewrite(xmg);
	//cout << xmg.num_gates() << " ";
	Refactor(xmg);
	//cout << xmg.num_gates() << " ";
	Balance(xmg);
	//cout << xmg.num_gates() << " ";
	Rewrite(xmg);
	//cout << xmg.num_gates() << " ";
	Rewrite(xmg, true);
	//cout << xmg.num_gates() << " ";
	Balance(xmg);
	//cout << xmg.num_gates() << " ";
	Refactor(xmg, true);
	//cout << xmg.num_gates() << " ";
	Rewrite(xmg, true);
	//cout << xmg.num_gates() << " ";
	Balance(xmg);
	//cout << xmg.num_gates() << "\n";
}