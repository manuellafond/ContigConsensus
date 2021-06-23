#pragma once
#include <algorithm>
#include <map>
#include <vector>

#include "Contig.h"
#include "CostMap.hpp"




using namespace std;


//Note: this is a struct, so everything is public
struct Match
{
  static const unsigned MATCHTYPE_PREFIX = 1;
  static const unsigned MATCHTYPE_SUFFIX = 2;
  static const unsigned MATCHTYPE_FULL_T_IN_S = 3;
  static const unsigned MATCHTYPE_FULL_S_IN_T = 4;

  struct Mcontig {
    Contig* contig;
    size_t start;
    bool isReverse;
  };

  array<Mcontig, 2> contigs;

  size_t length;
  unsigned score;

  Match() : score(0) {} // empty match
  Match(Contig *c1, Contig *c2, size_t start1, size_t start2, size_t length,
        bool isReverse1, bool isReverse2, CostMap &scores)
    : contigs{{{c1, start1, isReverse1},{c2,start2, isReverse2}}}, length(length), score(0) {


    for (int p = 0; p < length; ++p)
      score+= scores(c1->at(start1 + p, isReverse1),c2->at(start2 + p, isReverse2));
      
  }

  int GetType()
  {
    if (length == contigs[0].contig->size())
      return MATCHTYPE_FULL_S_IN_T;
    if (length == contigs[1].contig->size())
      return MATCHTYPE_FULL_T_IN_S;
    if (contigs[0].start == 0)
      return MATCHTYPE_PREFIX;
    if (contigs[1].start == 0)
      return MATCHTYPE_SUFFIX;
    throw std::runtime_error("This type of match is not supposed to happen.");
  }



	
  Contig * contig(unsigned t)  const
  {
    return contigs[t].contig;
  }

  Contig * contig(Contig * c) 
  {
    return c==contigs[0].contig ? contigs[1].contig : contigs[0].contig;
  }

  size_t start(unsigned t) const
  {
    return contigs[t].start;
  }

  size_t end(unsigned t) const
  {
    return start(t) + length-1;
  }

  bool is_reverse(unsigned t) const
  {
    return contigs[t].isReverse;
  }

  bool is_full(unsigned t)  const
  {
    return length == (contigs[t].contig->size());
  }

  size_t projected_start(unsigned t) const
  {
    return is_reverse(t) ? (contig(t)->size() - 1) - end(t) : start(t);
  }

  size_t projected_end(unsigned t) const
  {
    return is_reverse(t) ?  (contig(t)->size() - 1) - start(t) : end(t) ;
  }

  bool intersect(Match * m){
    if(this->contig(0)==m->contig(0) && this->contig(1)==m->contig(1))
      return true;
    
    for(unsigned i =0; i <=1; i++){
      if (this->contig(i) == m->contig(i) &&
	  ( this->is_reverse(i)!= m->is_reverse(i) ||
	    (this->start(i) >= m->start(i) && this->start(i) <= m->end(i)) ||
	    (this->end(i) >= m->start(i) && this->end(i) <= m->end(i)))) {
	return true;
      }
    }
    return false;
  }

  bool intersect(vector<Match*> v){
    return any_of(v.begin(), v.end(),
		   [&](Match *m) {return this->intersect(m);});
  }

  bool contains(Contig *c) // for debug
  {
    return contig(0)==c || contig(1)==c;
  }
  



	
};

ostream &operator<<(ostream& os, const Match& m)
{
  auto f = max((size_t)22,max(m.start(1)+m.contig(0)->size(),m.start(0)+m.contig(1)->size()));
  std::cout << endl << string(f,'-')<<endl;
  cout << "Len: " <<m.length<< " - Score: "<< m.score <<"\n";
  std::cout << string(f,'=')<<endl;

  for(int c : {0,1}){
    if(m.is_reverse(c)) cout << "R| ";
    else cout << " | ";
    std::cout << string(m.start((c+1)%2),' ');
    for (int i = 0; i < m.contig(c)->size(); i++) {
      if (i >= m.start(c) && i < m.start(c) + m.length) {
	if(m.contig(c)->at(i)==m.contig((c+1)%2)->at(i+m.start((c+1)%2)-m.start(c)))
	  cout << "\033[1;32m";
        else cout << "\033[1;31m";
      }

      cout << (char)m.contig(c)->at(i, m.is_reverse(c));
      cout << "\033[0m";
    }
    cout << endl;
  }
  std::cout << string(f,'=')<<endl;
  return os;
}



typedef map<Contig*,map<Contig*,Match*>> MM_map;
 
