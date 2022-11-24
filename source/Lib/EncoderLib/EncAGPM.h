
/** \file     EncAGPM.h
\brief
*/

#ifndef __ENCAGPM__
#define __ENCAGPM__

// Include files
#include "CommonLib/UnitTools.h"
#include "EncCu.h"
//! \ingroup EncoderLib
//! \{

 #if AGPM_ENABLED

class EncLib;
class HLSWriter;
class EncSlice;

// ====================================================================================================================
// Class definition
// ====================================================================================================================

/// CU encoder class
/*
struct AGPMMergeCombo
{
	int    splitDir;
	int    mergeIdx0;
	int    mergeIdx1;
	double cost; //RD Cost?
	AGPMMergeCombo() : splitDir(), mergeIdx0(-1), mergeIdx1(-1), cost(0.0) {};
	AGPMMergeCombo(int _splitDir, int _mergeIdx0, int _mergeIdx1, double _cost)
		: splitDir(_splitDir), mergeIdx0(_mergeIdx0), mergeIdx1(_mergeIdx1), cost(_cost) {};
};

struct AGPMMotionInfo
{
	uint8_t m_splitDir; //���ַ���
	uint8_t m_candIdx0; //��һ�������ĺ�ѡ���� 6��
	uint8_t m_candIdx1; //�ڶ��������ĺ�ѡ���� 6-1=5��

	AGPMMotionInfo(uint8_t splitDir, uint8_t candIdx0, uint8_t candIdx1) : m_splitDir(splitDir), m_candIdx0(candIdx0), m_candIdx1(candIdx1) { }
	AGPMMotionInfo() { m_splitDir = m_candIdx0 = m_candIdx1 = 0; }
};

struct AdaptivePartitionMotionInfo
{
	uint8_t m_candIdx0; //��һ�������ĺ�ѡ���� 6��
	uint8_t m_candIdx1; //�ڶ��������ĺ�ѡ���� 6-1=5��

	AdaptivePartitionMotionInfo(uint8_t candIdx0, uint8_t candIdx1) :m_candIdx0(candIdx0), m_candIdx1(candIdx1) { }
	AdaptivePartitionMotionInfo() {m_candIdx0 = m_candIdx1 = 0; }
};


struct BetterThanComboCost
{
	inline bool operator()(const AGPMMergeCombo &first, const AGPMMergeCombo &second) { return (first.cost < second.cost); }
};

class AGPMComboCostList
{
public:
	AGPMComboCostList() {};
	~AGPMComboCostList() {};
	std::vector<AGPMMergeCombo> list;

	void sortByCost() { std::stable_sort(list.begin(), list.end(), BetterThanComboCost()); };
};
*/
  #endif
#endif   // 
