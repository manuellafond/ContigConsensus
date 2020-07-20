#include <iostream>

#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include "Contig.h"
#include "Match.h"

using namespace std;







/**
* Filters matches to return only those that are type 1 (according to paper, prefix matches on non-reversed, 
* suffix matches on reversed, and full matches)
**/
vector<Match*> GetType1Matches(vector<Match*>& matches)
{
	vector<Match*> retval;

	for (int i = 0; i < matches.size(); ++i)
	{
		if ( (matches[i]->GetType() == Match::MATCHTYPE_PREFIX && !matches[i]->isReverse1) ||
			 (matches[i]->GetType() == Match::MATCHTYPE_SUFFIX && matches[i]->isReverse1)  ||
			 (matches[i]->GetType() == Match::MATCHTYPE_FULL_S_IN_T) ||
			 (matches[i]->GetType() == Match::MATCHTYPE_FULL_T_IN_S)
		   )
		{
			retval.push_back(matches[i]);
		}
	}

	return retval;
}




/**
* Filters matches to return only those that are type 2 (according to paper, prefix matches on reversed,
* suffix matches on non-reversed, and full matches)
**/
vector<Match*> GetType2MatchesIndices(vector<Match*>& matches)
{
	vector<Match*> retval;

	for (int i = 0; i < matches.size(); ++i)
	{
		if ((matches[i]->GetType() == Match::MATCHTYPE_PREFIX && matches[i]->isReverse1) ||
			(matches[i]->GetType() == Match::MATCHTYPE_SUFFIX && !matches[i]->isReverse1) ||
			(matches[i]->GetType() == Match::MATCHTYPE_FULL_S_IN_T) ||
			(matches[i]->GetType() == Match::MATCHTYPE_FULL_T_IN_S)
			)
		{
			retval.push_back(matches[i]);
		}
	}

	return retval;
}







//the set is passed as a copy because we will empty it
vector<Match*> GetBEO(set<Match*> allMatches)
{
	//NOT DONE YET
	vector<Match*> beoMatches;

	map<Contig*, set<Match*>> matchesPerContig;
	for (Match* m : allMatches)
	{
		matchesPerContig[m->contig1].insert(m);
		matchesPerContig[m->contig2].insert(m);
	}

	
	while (!allMatches.empty())
	{
		//find match with largest contig
		Match* largestMatch = *(allMatches.begin());
		int maxlen = max(largestMatch->contig1->size(), largestMatch->contig2->size());
		for (Match* m : allMatches)
		{
			int len = max(m->contig1->size(), m->contig2->size());
			if (len > maxlen)
			{
				largestMatch = m;
				maxlen = len;
			}
		}

		//find early ender inside largestMatch
		Contig* maxcontig = largestMatch->contig1;
		if (largestMatch->contig2->size() > maxcontig->size())
			maxcontig = largestMatch->contig2;

		Contig* earlyender = nullptr;
		int earlyenderpos = 99999;
		Contig* latestarter = nullptr;
		int latestarterpos = -1;
		for (Match* m : matchesPerContig[maxcontig])
		{
			int end, start;
			if (maxcontig == m->contig1)
			{
				start = m->start1;
				end = m->start1 + m->length - 1;
			}
			else
			{
				start = m->start2;
				end = m->start2 + m->length - 1;
			}


		}
	}


	return beoMatches;
}








int main()
{
	//for testing purposes
	vector<Contig*> S;
	S.push_back(new Contig({0,0,0,0,0,0,1,1,1,1}));

	vector<Contig*> T;
	T.push_back(new Contig({1, 0, 1, 0}));
	T.push_back(new Contig({ 0, 0, 0 }));


	//scores are just 0 if same char, 1 if different char
	//TODO: should actually be called costs, not scores
	map<int, map<int, int>> scores;

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			int score = 0;
			if (i == j)
				score = 0;
			scores[i][j] = 1;
		}
	}

	vector<Match*> matches;

	for (Contig* C : S)
	{
		for (Contig* D : T)
		{
			int minlen = min(C->size(), D->size());
			int maxlen = max(C->size(), D->size());

			//this loop computes suffix-prefix matches and prefix-suffix matches
			for (int i = 0; i < minlen; ++i)
			{
				//compute all suffix-prefix matches of length i + 1
				matches.push_back( new Match(C, D, C->size() - i - 1, 0, i + 1, false, false, scores) );
				matches.push_back(new Match(C, D, C->size() - i - 1, 0, i + 1, true, false, scores) );
				matches.push_back(new Match(C, D, C->size() - i - 1, 0, i + 1, false, true, scores));
				matches.push_back(new Match(C, D, C->size() - i - 1, 0, i + 1, true, true, scores));
				

				//compute all prefix-suffix matches of length i + 1
				matches.push_back(new Match(C, D, 1, D->size() - i - 1, i + 1, false, false, scores));
				matches.push_back(new Match(C, D, 1, D->size() - i - 1, i + 1, true, false, scores));
				matches.push_back(new Match(C, D, 1, D->size() - i - 1, i + 1, false, true, scores));
				matches.push_back(new Match(C, D, 1, D->size() - i - 1, i + 1, true, true, scores));
				
			}

			//compute all full matches
			if (C->size() <= D->size())
			{
				//here i has all starting pos in D
				for (int i = 0; i <= D->size() - C->size(); i++)
				{
					matches.push_back(new Match(C, D, 0, i, C->size(), false, false, scores) );
					matches.push_back(new Match(C, D, 0, i, C->size(), true, false, scores) );
					matches.push_back(new Match(C, D, 0, i, C->size(), false, true, scores) );
					matches.push_back(new Match(C, D, 0, i, C->size(), true, true, scores) );
				}
			}
			else
			{
				//maybe we could do it in one loop, but clarity might be lost
				//here i has all starting pos in D
				for (int i = 0; i <= C->size() - D->size(); i++)
				{
					matches.push_back(new Match(C, D, i, 0, D->size(), false, false, scores) );
					matches.push_back(new Match(C, D, i, 0, D->size(), true, false, scores) );
					matches.push_back(new Match(C, D, i, 0, D->size(), false, true, scores) );
					matches.push_back(new Match(C, D, i, 0, D->size(), true, true, scores) );
				}
			}
		}
	}



	//for debugging
	for (Match* m : matches)
	{
		cout << "Matching ";
		for (int i = 0; i < m->contig1->size(); ++i)
			cout << m->contig1->at(i);
		cout << "   vs   ";
		for (int i = 0; i < m->contig2->size(); ++i)
			cout << m->contig2->at(i);
		cout << endl;
		cout << "s1=" << m->start1 << " s2=" << m->start2 << " r1=" << m->isReverse1 
			 << " r2=" << m->isReverse2 << " l="<<m->length<<" s="<<m->score<<endl;
	}




	for (Match* m : matches)
	{
		delete m;
	}


	for (Contig* C : S)
	{
		delete C;
	}

	for (Contig* C : T)
	{
		delete C;
	}
}

