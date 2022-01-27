#include <iostream>
#include <stdexcept>
#include <string>
#include <time.h>
#include <unistd.h>
#include <utility>

#include "Contig.h"
#include "Hungarian_algo.hpp"
#include "Match.h"
#include "Parser.hpp"
#include "MatchesPerContig.hpp"
#include "Nucleotide.h"

using namespace std;

unsigned greedy_fill(vector<const Match*>& selected_matches, vector<Match>& all_matches){
  unsigned added_score=0;

  sort(all_matches.begin(),all_matches.end(),[](const Match &m1, const Match& m2){
    return m1.score>m2.score;
   });

  for(Match & m : all_matches){
    if (!m.intersect(selected_matches)) {
      selected_matches.push_back(&m);
      added_score += m.score;
    }
  }
  return added_score;
}




void update_match(Match &m,map<Contig*,tuple<Contig*,size_t,bool>>& m_contig, MatchMatrix& matches, unsigned other_set_id, unsigned new_set_id)
{

  size_t shift_c1=0, shift_c2=0;
  bool is_c1_reversed=false , is_c2_reversed=false;
  Contig * c1 = m.contig((unsigned)0);
  
  while(m_contig.find(c1)!=m_contig.end()){
    auto &r = m_contig[c1];
    c1=get<0>(r);
    shift_c1+=get<1>(r);
    is_c1_reversed=(is_c1_reversed!=get<2>(r));
  }
  unsigned start_c1= is_c1_reversed ? shift_c1 + m.relative_start((unsigned)0) : shift_c1 + m.start((unsigned)0);
    
  Contig * c2 = m.contig(1);
  while(m_contig.find(c2)!=m_contig.end()){
    auto &r = m_contig[c2];
    c2= get<0>(r);
    shift_c2+=get<1>(r);
    is_c2_reversed=(is_c2_reversed!=get<2>(r));
  }

  unsigned start_c2= is_c2_reversed ? shift_c2 + m.relative_start(1) : shift_c2 + m.start(1);

  get<0>(matches[other_set_id][new_set_id]).emplace_back(Match(c1,c2,start_c1,is_c1_reversed!=m.is_reverse((unsigned)0),
							       start_c2,is_c2_reversed!=m.is_reverse(1)));


}



vector<const Match*> merge_full_matches(vector<const Match *> &selected_matches,
                        AssemblySet &assembly_sets, MatchMatrix &matches,
			vector<string> &trash, map<Contig*,tuple<Contig*,size_t,bool>> &m_contig,
                        unsigned new_set_id, unsigned set_id1, unsigned set_id2)
{
  auto &new_set =assembly_sets[new_set_id];

  vector<const Match*> other_matches;
  for(auto &m : selected_matches){
    if(!m->is_full((unsigned) 0) && !m->is_full(1)){
      other_matches.push_back(m);
      continue;
    }
    
    Contig * c1 = m->contig((unsigned)0);
    unsigned shift_c1=0;
    while(m_contig.find(c1)!=m_contig.end())
      c1=get<0>(m_contig[c1]);

    Contig * c2 = m->contig((unsigned)1);
    unsigned shift_c2=0;
    while(m_contig.find(c2)!=m_contig.end())
      c2=get<0>(m_contig[c2]);


    unsigned overlap;
    bool is_c1_reversed, is_c2_reversed;
    if(!m->is_full(1)){
      overlap = m->is_reverse(1) ?  m->contig(1)->size()-m->relative_start(1) : m->contig(1)->size()-m->start(1);
      swap(c1,c2);
      is_c2_reversed=m->is_reverse((unsigned)0);
      is_c1_reversed=m->is_reverse(1);
    } else {
      overlap = m->is_reverse((unsigned)0) ?  m->contig((unsigned)0)->size()-m->relative_start((unsigned)0) : m->contig((unsigned)0)->size()-m->start((unsigned)0);

      is_c1_reversed=m->is_reverse((unsigned)0);
      is_c2_reversed=m->is_reverse(1);
    } 
  

      Contig * s = new_set.emplace(move(make_unique<Contig>(move(*c1),move(*c2),
							    is_c1_reversed,
							    is_c2_reversed,
							    overlap,
							    new_set_id))).first->get();
      

      m_contig[c1]=make_tuple(s,0,is_c1_reversed);
      m_contig[c2]=make_tuple(s, overlap, is_c2_reversed);

      if(c1->get_set_id()==new_set_id)
	trash.push_back(c1->getName());
      if(c2->get_set_id()==new_set_id)
	trash.push_back(c2->getName());
      
  }

  return other_matches;

  
}


void merge_other_matches(vector<const Match *> &selected_matches, vector<string> &trash,
			 AssemblySet &assembly_sets, MatchMatrix& matches,  map<Contig*,tuple<Contig*,size_t, bool>>& m_contig,
			 unsigned new_set_id,unsigned set_id1 = 1,
                 unsigned set_id2 = 2)
{
  auto &new_set =assembly_sets[new_set_id];

  for(auto &m : selected_matches){    
    Contig * c1 = m->contig((unsigned)0);


    unsigned shift_c1=0;
    bool is_c1_reversed=m->is_reverse((unsigned)0);
    while(m_contig.find(c1)!=m_contig.end()){
      auto &tmp = m_contig[c1];
      c1=get<0>(tmp);
      shift_c1+=get<1>(tmp);
      is_c1_reversed=(is_c1_reversed!=get<2>(tmp));
    }

    Contig * c2 = m->contig(1);
    unsigned shift_c2=0;
    bool is_c2_reversed=m->is_reverse(1);
    while(m_contig.find(c2)!=m_contig.end()){
      auto &tmp= m_contig[c2];
      c2= get<0>(tmp);
      shift_c2+=get<1>(tmp);
      is_c2_reversed=(is_c2_reversed!=get<2>(tmp));
    }


    //size_t old_size_c1 =c1->size(), old_size_c2=c2->size();
    unsigned relative_start_c1=is_c1_reversed ? c1->size()-1-m->projected_end((unsigned)0) : m->start((unsigned)0);
    unsigned relative_start_c2=is_c2_reversed ? c2->size()-1-m->projected_end(1) : m->start(1);

    if(relative_start_c1>relative_start_c2) {
      Contig * s = new_set.emplace(move(make_unique<Contig>(move(*c1),move(*c2),
							    is_c1_reversed,
							    is_c2_reversed,
							    m->length(1),new_set_id))).first->get();

      m_contig[c1]=make_tuple(s, 0, is_c1_reversed);
      m_contig[c2]=make_tuple(s, s->size()-c2->size(), is_c2_reversed);
      
      // cout << c1->getName() << "-> (" << s->getName() << ",0," <<  is_c1_reversed << ")" << endl;

      // cout << c2->getName() << "-> (" << s->getName() << ","<< s->size()-old_size_c2 << "," << is_c2_reversed << ")" << endl;
    
    } else {
      Contig * s = new_set.emplace(move(make_unique<Contig>(move(*c2),move(*c1),
							    is_c2_reversed,
							    is_c1_reversed,
							    m->length(1),
							    new_set_id))).first->get();


      // cout << c2->getName() << "-> (" << s->getName() << ",0," << is_c2_reversed << ")" << endl;

      // cout << c1->getName() << "-> (" << s->getName() << ","<<  s->size()-old_size_c1 << "," << is_c1_reversed << ")" << endl;

      m_contig[c2]=make_tuple(s, 0, is_c2_reversed);
      m_contig[c1]=make_tuple(s, s->size()-c1->size(), is_c1_reversed);
      
    }
    
    if(c1->get_set_id()==new_set_id)
      trash.push_back(c1->getName());
    if(c2->get_set_id()==new_set_id)
      trash.push_back(c2->getName());
  }

}




void merge_match(AssemblySet &assembly_sets, MatchMatrix& matches, unsigned new_set_id,unsigned set_id1,
                 unsigned set_id2)
{
  auto &selected_matches = get<2>(matches[set_id1][set_id2]);

  vector<string> trash;
  map<Contig *, tuple<Contig *, size_t, bool>> m_contig;
  auto v=  merge_full_matches(selected_matches, assembly_sets, matches, trash, m_contig, new_set_id, set_id1, set_id2);
  merge_other_matches(v, trash, assembly_sets, matches, m_contig, new_set_id, set_id1, set_id2);


  // update other matches
  auto &new_set =assembly_sets[new_set_id];
  for(auto &a : matches){
    if(a.first==set_id1){
      for(auto &b : a.second){
	if (b.first == set_id2)
          continue;
	for(auto &m : get<0>(b.second)){
	  update_match(m, m_contig, matches, b.first, new_set_id);
	}
      }
    }
    else if(a.first<set_id1){
      for(auto &m : get<0>(a.second[set_id1])){
	update_match(m, m_contig, matches, a.first, new_set_id);
      }
    }
    if(a.first==set_id2){
      for(auto &b : a.second){
        if (b.first == set_id1)
          continue;
	for(auto &m : get<0>(b.second)){
	  update_match(m, m_contig, matches, b.first, new_set_id);}
      }
    }
    else if(a.first<set_id2){
      for(auto &m : get<0>(a.second[set_id1])){
	update_match(m, m_contig, matches, a.first, new_set_id);
      }
    }
  }


  for(auto it = new_set.begin(); it!=new_set.end();){
    if(m_contig.find(it->get())!=m_contig.end()){
      // remove intermediate contigs
      it=new_set.erase(it);
    } else ++it;
    
  }

  assembly_sets.erase(set_id1);
  assembly_sets.erase(set_id2);
  

  matches.erase(set_id1);
  matches.erase(set_id2);

  for(auto &a : matches){
    a.second.erase(set_id1);
    a.second.erase(set_id2);
  }
}

void merge_algorithm(AssemblySet &assembly_sets, MatchMatrix &matches, map<string,unsigned> &ids){

  
  for(size_t i=0; i<ids.size();++i)
    for(size_t j=i+1;j<ids.size();++j){
      auto &m = matches[i][j];
      get<1>(m) =  greedy_fill(get<2>(m),get<0>(m));
    }

  unsigned new_id=ids.size();
  unsigned cpt_id=ids.size();
  for(unsigned cpt=0;cpt<cpt_id-1;++cpt){
  
    unsigned max_score=0;
    unsigned i=matches.begin()->first,
      j=matches.begin()->second.begin()->first;
    
    // select two sets with a maximum score
    for(auto &a : matches)
      for(auto &b : a.second) {
	if(max_score<= get<1>(b.second)){
	  max_score=get<1>(b.second);
	  i=a.first;
	  j=b.first;
	}
      }

    //    cout << "Merge " << i << " with " << j << endl;

    merge_match(assembly_sets, matches, new_id, i ,j);
    
    
      for(auto & a : matches){
	auto &m = matches[a.first][new_id];
	get<1>(m) =  greedy_fill(get<2>(m), get<0>(m));//algo(assembly_sets, get<0>(m));

      }
      new_id++;
  }
   
}
