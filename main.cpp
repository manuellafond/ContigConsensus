#include <iostream>

#include <algorithm>
#include <map>
#include <set>
#include <vector>

#include <unistd.h>


#include "Match.h"
#include "Parser.hpp"

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


void algo(AssemblySet &T, AssemblySet &S, CostMap &scores)
{
	//for testing purposes
	

	//scores are just 0 if same char, 1 if different char
	//TODO: should actually be called costs, not scores
	//CostMap scores;

	vector<Match*> matches;

	for (Contig& C : S)
	{
		for (Contig& D : T)
		{
			int minlen = min(C.size(), D.size());
			int maxlen = max(C.size(), D.size());

			//this loop computes suffix-prefix matches and prefix-suffix matches
			for (int i = 0; i < minlen; ++i)
			{
				//compute all suffix-prefix matches of length i + 1
				matches.push_back( new Match(&C, &D, C.size() - i - 1, 0, i + 1, false, false, scores) );
				matches.push_back(new Match(&C, &D, C.size() - i - 1, 0, i + 1, true, false, scores) );
				matches.push_back(new Match(&C, &D, C.size() - i - 1, 0, i + 1, false, true, scores));
				matches.push_back(new Match(&C, &D, C.size() - i - 1, 0, i + 1, true, true, scores));
				

				//compute all prefix-suffix matches of length i + 1
				matches.push_back(new Match(&C, &D, 1, D.size() - i - 1, i + 1, false, false, scores));
				matches.push_back(new Match(&C, &D, 1, D.size() - i - 1, i + 1, true, false, scores));
				matches.push_back(new Match(&C, &D, 1, D.size() - i - 1, i + 1, false, true, scores));
				matches.push_back(new Match(&C, &D, 1, D.size() - i - 1, i + 1, true, true, scores));
				
			}

			//compute all full matches
			if (C.size() <= D.size())
			{
				//here i has all starting pos in D
				for (int i = 0; i <= D.size() - C.size(); i++)
				{
					matches.push_back(new Match(&C, &D, 0, i, C.size(), false, false, scores) );
					matches.push_back(new Match(&C, &D, 0, i, C.size(), true, false, scores) );
					matches.push_back(new Match(&C, &D, 0, i, C.size(), false, true, scores) );
					matches.push_back(new Match(&C, &D, 0, i, C.size(), true, true, scores) );
				}
			}
			else
			{
				//maybe we could do it in one loop, but clarity might be lost
				//here i has all starting pos in D
				for (int i = 0; i <= C.size() - D.size(); i++)
				{
					matches.push_back(new Match(&C, &D, i, 0, D.size(), false, false, scores) );
					matches.push_back(new Match(&C, &D, i, 0, D.size(), true, false, scores) );
					matches.push_back(new Match(&C, &D, i, 0, D.size(), false, true, scores) );
					matches.push_back(new Match(&C, &D, i, 0, D.size(), true, true, scores) );
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

}



array<string, 4> treatProgrammeEntry(int argc, char * argv[])
{
	int opt;
	array<string, 4> optionsValues = {"-1","-1","-1","-1"};

	while((opt = getopt(argc, argv, ":hs:t:c:o:")) !=-1)
	{
		switch (opt)
		{
		case 't': 
			optionsValues[0]=optarg;
			break;
		case 's':
			optionsValues[1]=optarg;
			break;
		
		case 'c': 
			optionsValues[2]=optarg;
			break;
		
		case 'h':
		default:
			cout << "Available options:\n";
			cout << "-t fasta file containing the T set\n";
			cout << "-s fasta file containing the S set\n";
			cout << "-o output fasta file name (not implemented)\n";
			cout << "-c cost map file (optionnal)\n";
			break;
		}
	}
	if(optionsValues[0]=="-1")
	{
		cout << "Error, \"T\" file not provided\n";
		exit(EXIT_FAILURE);
	}
	if(optionsValues[1]=="-1")
	{
		cout << "Error, \"S\" file not provided\n";
		exit(EXIT_FAILURE);
	}


	return optionsValues;
}



// test: ./ContigConsensus -t test/t1.fasta -s test/s1.fasta
int main(int argc, char *argv[])
{
	auto options = treatProgrammeEntry(argc, argv);
	CostMap c(options[2].c_str());
	
	cout << "Reading T input file: " << options[0] << endl;
	AssemblySet T = parseFile(options[0].c_str());
	cout << "Done\n\n";
	cout << "Reading S input file: " << options[1] << endl;
	AssemblySet S = parseFile(options[1].c_str());
	cout << "Done\n\n";

	algo(T, S,c);
}

