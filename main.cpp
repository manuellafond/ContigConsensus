#include <iostream>
#include <stdexcept>
#include <string>
#include <time.h>
#include <unistd.h>

#include "Hungarian_algo.hpp"
#include "Match.h"
#include "Parser.hpp"
#include "MatchesPerContig.hpp"


using namespace std;

queue<const Match*> construct_BEO(Match::SortedLengthSet& matches)
{

  queue<const Match*> beoMatches;
  unsigned nb_non_full_contigs=0;
  map<const Contig*,MatchesPerContig> mpc;
  set<const Match*> non_full_matches;

  for (const Match *m : matches){
    mpc.emplace(m->contig((unsigned)0), MatchesPerContig(0));
    mpc.emplace(m->contig(1), MatchesPerContig(1));

    nb_non_full_contigs+=mpc.at(m->contig((unsigned)0)).add(m);
    nb_non_full_contigs+=mpc.at(m->contig(1)).add(m);
  }
  while (!matches.empty()){
    if(nb_non_full_contigs==0){
      for (const Match *m : matches) 
        beoMatches.push(m);
      
      return beoMatches;
    }

		
    auto p = *find_if(mpc.begin(),mpc.end(),[](const pair<const Contig *, MatchesPerContig> &p){
      return !p.second.only_full();
    });

            
    const Contig * cI = p.first;
    const Match* IJ = p.second.early();
    const Contig *cJ; // equivalent to w

    while(true){

      cJ = IJ->contig(cI);
	
      if (mpc.at(cJ).is_early_or_late(IJ)) 
	break;
      

      if(IJ->projected_end(cJ->get_set_id())==cJ->size()-1){
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
    matches.erase(IJ);
  }
  return beoMatches;
}

vector<const Match*> local_ratio_rec(queue<const Match*> &beo, vector<pair<const Match*,int>> &sigma, unsigned& score)
{
  if(beo.empty())
    return vector<const Match*>();
  
  const  Match * F = beo.front();
  beo.pop();

  unsigned score_f = F->score;
  for(auto &p : sigma){
    if(F->intersect(p.first)) // what if p.second<0?
      score_f-= p.second;
  }

  if(score_f>0)
    sigma.push_back(make_pair(F, score_f));
  
  vector<const Match*> M = local_ratio_rec(beo,sigma,score);
  
  if (score_f>0 && !F->intersect(M)) {
    M.push_back(F);
    score+=F->score;
  }

  return M;
}


vector<const Match*> local_ratio(queue<const Match*> &beo, unsigned& score)
{
  vector<pair<const Match*,int>> sigma;
  return local_ratio_rec(beo,sigma,score);
}



unsigned greedy_fill(vector<const Match*>& selected_matches, vector<Match>& all_matches){
  unsigned added_score=0;

  for(Match & m : all_matches){
    if (!m.intersect(selected_matches)) {
      selected_matches.push_back(&m);
      added_score += m.score;
    }
  }
  return added_score;
}


unsigned random_greedy_fill(vector<const Match*>& selected_matches, vector<Match>& all_matches, unsigned probability=90){
  unsigned added_score=0;
  srand (time(NULL));
  for(Match & m : all_matches){
    if((rand()%100)>probability)
      continue;
    if (!m.intersect(selected_matches)) {
      selected_matches.push_back(&m);
      added_score += m.score;
    }
  }
  return added_score;
}

vector<const Match*> greedy_maximum_matching(const Match::MM_map& matches, unsigned &score)
{
  score=0;
  vector<const Match*> selected_matches;
  for(auto p1 : matches)
    for(auto p2 : p1.second)
      selected_matches.push_back(p2.second);

  sort(selected_matches.begin(),selected_matches.end(),[](const Match *m1, const Match* m2){
    return m1->score>m2->score;
  });

  vector<const Match*> rtrn;
  for(const Match * m : selected_matches){
    if (!m->intersect(rtrn)) {
      rtrn.push_back(m);
      score += m->score;
    }
  }
  return rtrn;
}




vector<const Match*> algo(AssemblySet &assembly_sets, vector<Match> &all_matches, unsigned id_set1=1, unsigned id_set2=2)
{
  Match::MM_map max_matches;
  unsigned maximum_value=0;


  cout << "Compute Greedy" << endl;
  vector<const Match*> selected_matches_best;
  cout << "Sort:\n";
  sort(all_matches.begin(),all_matches.end(),[](const Match &m1, const Match& m2){
    return m1.score>m2.score;
  });
  
  unsigned score_max = greedy_fill(selected_matches_best, all_matches);
  cout << "Score greedy: " <<score_max <<endl;

  unsigned score_tmp;
  vector<const Match*> selected_matches_tmp;
  
  for(unsigned i=0;i<2;++i){
    cout << "Filter matches (type " << i <<")\n";
    auto m = Match::match_filter(i, all_matches);
    cout << "Done!\n";

    cout << "Construct BEO\n";
    auto beo = construct_BEO(m);
    cout << "Done!\n";

    score_tmp=0;
    selected_matches_tmp = local_ratio(beo,score_tmp);

    cout << "Score: " << score_tmp <<endl;
    score_tmp+= greedy_fill(selected_matches_tmp, all_matches);
    cout << "Score after greedy fill: " << score_tmp <<endl;

    if(score_tmp>score_max){
      score_max=score_tmp;
      selected_matches_best=move(selected_matches_tmp);
    }
  }
  
  
  for(Match &m : all_matches){
    try{
      if(max_matches.at(m.contig((unsigned)0)).at(m.contig(1))->score < m.score) {
	max_matches[m.contig((unsigned)0)][m.contig(1)] = &m;
	maximum_value = max(maximum_value, m.score);
      }
    }
    catch(out_of_range) {
      max_matches[m.contig((unsigned)0)][m.contig(1)] = &m;
      maximum_value = max(maximum_value, m.score);
    }
  }
  
  score_tmp=0;
  
  selected_matches_tmp = greedy_maximum_matching(max_matches, score_tmp);  //hungarian_algorithm(S, T, max_matches, maximum_value,score_m3);
  cout << "Hungarian algorithm done\n";

  cout << "Score maximum matching: " << score_tmp <<endl;
  score_tmp+= greedy_fill(selected_matches_tmp, all_matches);
  cout << "\nScore 3 after greedy fill: " << score_tmp <<endl;

  if(score_tmp>score_max){
    score_max=score_tmp;
    selected_matches_best=move(selected_matches_tmp);
  }


  return selected_matches_best;
}

void merge_match(vector<const Match *> &selected_matches,
                 AssemblySet &assembly_sets, MatchMatrix& matches, unsigned new_set_id,unsigned set_id1 = 1,
                 unsigned set_id2 = 2)
{
  vector<string> trash;
  // for debug;
  size_t n = 0;

  auto &new_set =assembly_sets[new_set_id];

  map<Contig*,Contig*> m_contig;

  for(auto &m : selected_matches){
    cout << "merge " << m->contig((unsigned)0)->getName() << " with " << m->contig(1)->getName() << endl;
    Contig * c1 = m->contig((unsigned)0);
    while(m_contig.find(c1)!=m_contig.end())
      c1=m_contig[c1];
    
    Contig * c2 = m->contig(1);
    while(m_contig.find(c2)!=m_contig.end())
      c2=m_contig[c2];

    cout << "merge " << c1->getName() << " with " << c2->getName() << endl;
    Contig * s = c1->size()==m->length((unsigned)0) ? //otherwise the shift calculation is wrong
      new_set.emplace(move(make_unique<Contig>(move(*c2),move(*c1),
					       m->is_reverse((unsigned )0)||m->is_reverse(1),
					       m->length(1),to_string(n),
					       new_set_id))).first->get() :
      new_set.emplace(move(make_unique<Contig>(move(*c1),move(*c2),
					       m->is_reverse((unsigned )0)||m->is_reverse(1),
					       m->length(1),to_string(n),new_set_id))).first->get();
    n++;
    m_contig[c1]=s;
    m_contig[c2]=s;
    cout << s->size() << endl;
    if(c1->get_set_id()==new_set_id)
      trash.push_back(c1->getName());
    if(c2->get_set_id()==new_set_id)
      trash.push_back(c2->getName());


    // // update other matches
    // vector<Contig*> prefix_or_full_c1, suffix_c1, suffix_or_full_c2, prefix_c2;
    


    // // update other matches
    // for(auto &a : matches){
    //   for(auto &b : a.second){
    // 	if(b.first==new_set_id) continue;
    // 	unsigned other_set_id;
    // 	if(a.first==set_id1 || a.first != set_id2)
    // 	  other_set_id=b.first;
    // 	else if  (b.first==set_id1 || b.first == set_id2)
    // 	  other_set_id=a.first;
    // 	else continue;
    // 	auto &v = b.second;
    // 	for(auto &m_c1 : v){
    // 	  if(!m_c1.contains(c1)) continue;
	  
    // 	  const Contig *c3= m_c1.contig(c1);
    // 	  for(auto &m_c2 :v){
    // 	    if(!m_c2.contains(c2) || !m_c2.contains(c3)) continue;
	    
    // 	  }
    // 	}

    //   }
    // }


    
  }

  for(auto it = new_set.begin(); it!=new_set.end();){
    if(m_contig.find(it->get())!=m_contig.end()){
      // remove intermediate contigs
      it=new_set.erase(it);
    } else ++it;
    
  }

  assembly_sets.erase(set_id1);
  assembly_sets.erase(set_id2);
  
  
}



array<string, 2> treatProgrammeEntry(int argc, char * argv[])
{
  int opt;
  array<string, 2> optionsValues = {"-1","output.cons"};

  while((opt = getopt(argc, argv, ":hi:o:")) !=-1) {
      switch (opt)
	{
	case 'i': 
	  optionsValues[0]=optarg;
	  break;
       	case 'o': 
	  optionsValues[1]=optarg;
	  break;
		
	case 'h':
	default:
	  cout << "Available options:\n";
	  cout << "-o output .cons file name (default output.cons)\n";
	  cout << "-i input file name \n";
	  cout << "The input file must be generated by blastn using the command: \n";
	  cout << "\t blastn -task megablast -query file1.fasta -subject file2.fasta -out inputFileName -outfmt \"6 score qseqid qstart qend qlen sseqid sstart send slen\"";
	  cout << endl;
	  exit(EXIT_SUCCESS);
	  break;
	}
  }
  if(optionsValues[0]=="-1") {
      cout << "Error, input file not provided\n";
      exit(EXIT_FAILURE);
  }

  return optionsValues;
}


int main(int argc, char *argv[])
{
  AssemblySet assembly_sets;
  MatchMatrix matches;

  auto options = treatProgrammeEntry(argc, argv);
  parseMatches(options[0].c_str(), assembly_sets, matches);
  
  vector<const Match*> result =  algo(assembly_sets, matches[1][2]);

  unsigned new_id_set=3;
  merge_match(result, assembly_sets,matches,new_id_set);


  exit(EXIT_SUCCESS);
  output(options[1].c_str(),result);
}

