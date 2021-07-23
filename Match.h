#pragma once
#include <algorithm>
#include <map>
#include <set>
#include <vector>

#include "Contig.h"
#include "CostMap.hpp"




using namespace std;


//Note: this is a struct, so everything is public
struct Match
{
  static bool shorter(const Match *m1,const Match *m2)
  {
    return m1->length(0)+m1->length(1) < m2->length(0)+m2->length(1);
  }

  typedef set<const Match*,decltype(&shorter)> SortedLengthSet;
  typedef map<const Contig*,map<const Contig*,Match*>> MM_map;


  /**
   * Filters matches to return only those that are type 1 or 2 
   **/
  static SortedLengthSet match_filter(unsigned t, vector<Match> & matches)
  {
    SortedLengthSet retval(&Match::shorter);

    for (int i = 0; i < matches.size(); ++i) 
      if (matches[i].is_of_type(t))
	retval.insert(&matches[i]);
    return retval;
  }

  struct Mcontig {
    const Contig* contig;
    size_t start, end;
  };

  array<Mcontig, 2> contigs;

  //size_t length;
  unsigned score;

  Match() : score(0) {} // empty match
  Match(const Contig *c1, const Contig *c2, size_t start1, size_t end1, size_t start2, size_t end2, unsigned score)
    : contigs{{{c1, start1, end1},{c2,start2, end2}}}, score(score) {}

  bool is_of_type(unsigned t) const
  {
    return projected_start(t)==0
      || (projected_start((t+1)%2)==0 && projected_end((t+1)%2)==contig((t+1)%2)->size());
  }

  const Contig * contig(unsigned t)  const
  {
    return contigs[t].contig;
  }

  const Contig * contig(const Contig * c) const
  {
    return c==contigs[0].contig ? contigs[1].contig : contigs[0].contig;
  }

  size_t start(unsigned t) const
  {
    return contigs[t].start;
  }

  size_t end(unsigned t) const
  {
    return contigs[t].end;
  }

  size_t length(unsigned t) const
  {
    return is_reverse(t) ? start(t) - end(t) : end(t) - start(t);
  }

  bool is_reverse(unsigned t) const
  {
    return contigs[t].start > contigs[t].end;
  }

  bool is_full(unsigned t)  const
  {
    return length(t) == (contigs[t].contig->size());
  }

  size_t projected_start(unsigned t) const
  {
    // TODO is that just min(start,end)?
    return min(start(t),end(t));
  }

  size_t projected_end(unsigned t) const
  {
    return max(start(t),end(t));
  }

  bool intersect(const Match * m) const
  {
    if(this->contig((unsigned)0)==m->contig((unsigned)0) && this->contig(1)==m->contig(1))
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

  bool intersect(vector<const Match*> v) const
  {
    return any_of(v.begin(), v.end(),
		   [&](const Match *m) {return this->intersect(m);});
  }

  bool contains(const Contig *c) const // for debug
  {
    return contig((unsigned)0)==c || contig(1)==c;
  }

};

// ostream &operator<<(ostream& os, const Match& m)
// {
//   auto f = max((size_t)22,max(m.start(1)+m.contig(0)->size(),m.start(0)+m.contig(1)->size()));
//   std::cout << endl << string(f,'-')<<endl;
//   cout << /*"Len: " <<m.length<<*/ " - Score: "<< m.score <<"\n";
//   std::cout << string(f,'=')<<endl;

//   for(int c : {0,1}){
//     if(m.is_reverse(c)) cout << "R| ";
//     else cout << " | ";
//     std::cout << string(m.start((c+1)%2),' ');
//     for (int i = 0; i < m.contig(c)->size(); i++) {
//       if (i >= m.start(c) && i < m.start(c) + m.length) {
// 	if(m.contig(c)->at(i)==m.contig((c+1)%2)->at(i+m.start((c+1)%2)-m.start(c)))
// 	  cout << "\033[1;32m";
//         else cout << "\033[1;31m";
//       }

//       cout << (char)m.contig(c)->at(i, m.is_reverse(c));
//       cout << "\033[0m";
//     }
//     cout << endl;
//   }
//   std::cout << string(f,'=')<<endl;
//   return os;
// }






