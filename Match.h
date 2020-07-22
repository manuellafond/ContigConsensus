#pragma once

#include "Contig.h"
#include "CostMap.hpp"
#include <map>

using namespace std;


//Note: this is a struct, so everything is public
struct Match
{
	static const int MATCHTYPE_PREFIX = 0;
	static const int MATCHTYPE_SUFFIX = 1;
	static const int MATCHTYPE_FULL_T_IN_S = 2;
	static const int MATCHTYPE_FULL_S_IN_T = 3;

	Contig* contig1;
	Contig* contig2;

	int start1;
	int start2;
	int length;
	bool isReverse1;
	bool isReverse2;

	int score;

	Match()
	{
		contig1 = nullptr;
		contig2 = nullptr;

		score = 0;
	}

	Match(Contig* c1, Contig* c2, int start1, int start2, int length, bool isReverse1, bool isReverse2, CostMap& scores)
	{
		this->contig1 = c1;
		this->contig2 = c2;
		this->start1 = start1;
		this->start2 = start2;
		this->isReverse1 = isReverse1;
		this->isReverse2 = isReverse2;
		this->length = length;

		this->score = 0;
		for (int p = 0; p < length; ++p)
		{
			int char1 = c1->at(start1 + p, isReverse1);
			int char2 = c2->at(start2 + p, isReverse2);
			this->score += scores(char1,char2);
		}

	}




	int GetType()
	{
		if (length == contig1->size())
			return MATCHTYPE_FULL_S_IN_T;
		else if (length == contig2->size())
			return MATCHTYPE_FULL_T_IN_S;
		else if (start1 == 0)
			return MATCHTYPE_PREFIX;
		else if (start2 == 0)
			return MATCHTYPE_SUFFIX;
		else
		{
			throw std::runtime_error("This type of match is not supposed to happen.");
		}
	}

	


	int GetEndPos1()
	{
		return start1 + length - 1;
	}

	int GetEndPos2()
	{
		return start2 + length - 1;
	}


	int GetProjectedStart1()
	{
		if (!isReverse1)
			return start1;
		else
			return (contig1->size() - 1) - GetEndPos1();
	}

	int GetProjectedEnd1()
	{
		if (!isReverse1)
			return GetEndPos1();
		else
			return (contig1->size() - 1) - start1;
	}


	int GetProjectedStart2()
	{
		if (!isReverse2)
			return start2;
		else
			return (contig2->size() - 1) - GetEndPos2();
	}

	int GetProjectedEnd2()
	{
		if (!isReverse2)
			return GetEndPos2();
		else
			return (contig2->size() - 1) - start2;
	}


	

};

