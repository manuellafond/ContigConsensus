#include <iostream>

#include <unistd.h>

#include "Hungarian_algo.hpp"
#include "Parser.hpp"
#include "MatchesPerContig.hpp"


using namespace std;


/**
 * Filters matches to return only those that are type 1 (according to paper, prefix matches on non-reversed, 
 * suffix matches on reversed, and full matches)
 **/
set<Match*> GetType1Matches(vector<Match>& matches)
{
  set<Match*> retval;

  for (int i = 0; i < matches.size(); ++i)
    {
      if ( (matches[i].GetType() == Match::MATCHTYPE_PREFIX && !matches[i].is_reverse(0)) ||
	   (matches[i].GetType() == Match::MATCHTYPE_SUFFIX && matches[i].is_reverse(0))  ||
	   (matches[i].GetType() == Match::MATCHTYPE_FULL_S_IN_T) ||
	   (matches[i].GetType() == Match::MATCHTYPE_FULL_T_IN_S)
	   )
	{
	  retval.insert(&matches[i]);
	}
    }

  return retval;
	
}




/**
 * Filters matches to return only those that are type 2 (according to paper, prefix matches on reversed,
 * suffix matches on non-reversed, and full matches)
 **/
set<Match*> GetType2MatchesIndices(vector<Match>& matches)
{
  set<Match*> retval;

  for (int i = 0; i < matches.size(); ++i)
    {
      if ((matches[i].GetType() == Match::MATCHTYPE_PREFIX && matches[i].is_reverse(0)) ||
	  (matches[i].GetType() == Match::MATCHTYPE_SUFFIX && !matches[i].is_reverse(0)) ||
	  (matches[i].GetType() == Match::MATCHTYPE_FULL_S_IN_T) ||
	  (matches[i].GetType() == Match::MATCHTYPE_FULL_T_IN_S)
	  )
	{
	  retval.insert(&matches[i]);
	}
    }

  return retval;
}







//the set is passed as a copy because we will empty it
queue<Match*> GetBEO(set<Match*> allMatches)
{

  queue<Match*> beoMatches;
	
  int nb_non_full_contigs=0;
	

  map<Contig*,MatchesPerContig> mpc;
  set<Match*> non_full_matches;



  for (Match *m : allMatches){
      
    mpc.insert({m->contig(0), MatchesPerContig(0)});
    mpc.insert({m->contig(1), MatchesPerContig(1)});

    
    nb_non_full_contigs+=mpc.at(m->contig(0)).add(m);
    nb_non_full_contigs+=mpc.at(m->contig(1)).add(m);

  }
  

  while (!allMatches.empty()){
      
    if(nb_non_full_contigs==0){
      for(Match *m : allMatches)
	beoMatches.push(m);
      return beoMatches;
    }

		
    auto p = *find_if(mpc.begin(),mpc.end(),[](const pair<Contig *, MatchesPerContig> &p){
      return !p.second.only_full();
    });

            
    Contig * cI = p.first;
    Match* IJ = p.second.early();
    Contig *cJ; // equivalent to w

    while(true){
		
      cJ = IJ->contig(cI);
	
	
      if (mpc.at(cJ).is_early_or_late(IJ)) {
	break;
      }

      if(IJ->projected_end(cJ->t)==cJ->size()-1){
	cI=cJ; 
	IJ= mpc.at(cI).late();
      }
      else {
	cI=cJ;
	IJ= mpc.at(cI).early();
      }
    }
      
    beoMatches.push(IJ);

    nb_non_full_contigs-= mpc.at(cI).remove(IJ);
    nb_non_full_contigs-= mpc.at(cJ).remove(IJ);
    allMatches.erase(IJ);

      
  }
	

  return beoMatches;
}


vector<Match*> maxMatches(queue<Match*> &beo, unsigned& score)
{
  if(beo.empty())
    return vector<Match*>();
  
  Match * F = beo.front();
  beo.pop();
  vector<Match*> M = maxMatches(beo,score);

  if (!F->intersect(M)) {
    M.push_back(F);
    score+=F->score;
  }

  return M;
}

unsigned greedy_fill(vector<Match*>& selected_matches, vector<Match>& all_matches){
  unsigned added_score=0;
  for(Match & m : all_matches)
    if (!m.intersect(selected_matches)) {
      selected_matches.push_back(&m);
      added_score+=m.score;
    }
  return added_score;
}



void algo(AssemblySet &T, AssemblySet &S, CostMap &scores, unsigned filter=0)
{
  //for testing purposes
  
  //scores are just 0 if same char, 1 if different char
  //TODO: should actually be called costs, not scores
  //CostMap scores;
  
  vector<Match> all_matches;
  MM_map max_matches;
  unsigned maximum_value=0;

  for (Contig& C : S) {
      for (Contig& D : T) {
	
	  max_matches[&C][&D]=nullptr;
	  int minlen = min(C.size(), D.size());
	  int maxlen = max(C.size(), D.size());

	  for(bool C_bool = false; !C_bool; C_bool=!C_bool)
	    for(bool D_bool = false; !D_bool; D_bool=!D_bool){

	      //this loop computes suffix-prefix matches and prefix-suffix matches
	      for (int i = 0; i < minlen; ++i)
		{
		  //compute all suffix-prefix matches of length i + 1
		  Match m = Match(&C, &D, C.size() - i - 1, 0, i + 1, C_bool, D_bool, scores);

		  
		  if(m.score>filter) {
		    all_matches.push_back(move(m));

		    if (max_matches[&C][&D] == nullptr ||
			max_matches[&C][&D]->score < all_matches.back().score) {
		      max_matches[&C][&D] = &all_matches.back();
		      maximum_value=max(maximum_value,max_matches[&C][&D]->score);
		    }
                  }
				
				

		  //compute all prefix-suffix matches of length i + 1
		  m = Match(&C, &D, 0, D.size() - i - 1, i + 1, C_bool, D_bool, scores);
		  if(m.score>filter) {
		    all_matches.push_back(move(m));

		    if (max_matches[&C][&D] == nullptr ||
			max_matches[&C][&D]->score < all_matches.back().score) {
		      max_matches[&C][&D] = &all_matches.back();
		      maximum_value=max(maximum_value,max_matches[&C][&D]->score);
		    }
                  }

				
				
				
		}
	      //compute all full matches
	      if (C.size() <= D.size())
		{
		  //here i has all starting pos in D
		  for (int i = 0; i <= D.size() - C.size(); i++) {
		    Match m = Match(&C, &D, 0, i, C.size(), C_bool, D_bool, scores);
		    if(m.score>0)
		      all_matches.push_back(move(m));
		  }
		}
	      else
		{
		  //maybe we could do it in one loop, but clarity might be lost
		  //here i has all starting pos in D
		  for (int i = 0; i <= C.size() - D.size(); i++) {

		    Match m = Match(&C, &D, i, 0, D.size(), C_bool, D_bool, scores);
		    if (m.score > filter)
		      all_matches.push_back(move(m));
		  }
		}

	      
            }
        }
    }
 

  
 
  
 auto m1 = GetType1Matches(all_matches);
 auto r = GetBEO(m1);

 cout << "BEO 1 done \n";

 unsigned score_1=0;
 auto s1 = maxMatches(r,score_1);

 
 // for(Match * m : s1){
 //   cout << *m;
 // }

 cout << "Score 1: " << score_1 <<endl;
 score_1+= greedy_fill(s1, all_matches);
 cout << "Score 1 after greedy fill: " << score_1 <<endl;

 auto m2 = GetType2MatchesIndices(all_matches);
 auto r2 = GetBEO(m2);

 cout << "BEO 2 done \n";
 unsigned score_2=0; 
 auto s2 = maxMatches(r2,score_2);


 cout << "Score 2: " << score_2 <<endl;
 score_2+= greedy_fill(s2, all_matches);
 cout << "Score 2 after greedy fill: " << score_2 <<endl;

 if(score_1>score_2){
   for(Match * m : s1)
   cout << *m;
 } else {
    for(Match * m : s2)
      cout << *m;
 }

  // TODO correct hungarian algo
 /*
  vector<Match*> selected_suf_pref_matches = hungarian_algorithm(T, S, max_matches, maximum_value);
    cout << "Matches build" << "\n";

  for(Match * m : selected_suf_pref_matches) {
    cout << *m;
  } 
 */
  

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
  AssemblySet T = parseFile(options[0].c_str(),0);
  cout << "Done\n\n";
  cout << "Reading S input file: " << options[1] << endl;
  AssemblySet S = parseFile(options[1].c_str(),1);
  cout << "Done\n\n";

  algo(T, S,c);
}

