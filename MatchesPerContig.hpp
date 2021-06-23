#pragma once

#include <memory>
#include <set>
#include <functional>

#include "Match.h"


using namespace std;


class MatchesPerContig
{
  Contig *c;
  int t;

  bool comp_early_ender(const Match* m1,const Match*m2){
    return m1->projected_end(t) == m2->projected_end(t) ?
      m1<m2 : m1->projected_end(t) < m2->projected_end(t) ;
  }

  bool comp_late_starter(const Match* m1, const Match*m2){
    return m1->projected_start(t) == m2->projected_start(t) ?
      m1<m2 : m1->projected_start(t) > m2->projected_start(t);
  }
 
  set<Match*,std::function<bool(const Match*, const Match*)>> early_ender;
  set<Match*,std::function<bool(const Match*, const Match*)>> late_starter;
  set<Match*> full;

public:
  MatchesPerContig(int t) : t(t),
			    early_ender(std::bind(&MatchesPerContig::comp_early_ender,
						  this,                                               
						  std::placeholders::_1,
						  std::placeholders::_2)),
			    late_starter(std::bind(&MatchesPerContig::comp_late_starter,
						   this,                                               
						   std::placeholders::_1,
						   std::placeholders::_2))
  {
    
  }

  // return true the first time a non-full interval is inserted
  bool add(Match *m){


    if (!m->is_full(t)) {
      if(!early_ender.insert(m).second) {cout << "error\n"; exit(EXIT_SUCCESS);};
      if(!late_starter.insert(m).second) {cout << "error\n"; exit(EXIT_SUCCESS);};
      return early_ender.size() == 1;
    }
    full.insert(m);
    return false;
  }


  // return true the first time all non-full interval have been deleted
  bool remove(Match *m)
  {
    
    if(!this->contains(m)) {
      cout << "Error " << "\n";
      exit(EXIT_SUCCESS);

    }
    //cout << "Remove" << endl;
    if(!m->is_full(t)){
      early_ender.erase(m);
      late_starter.erase(m);
      return early_ender.empty();
    }
    
    full.erase(m);
    return false;
  }

  bool contains(Match * m) // for debug purpose
  {
    return early_ender.count(m) || late_starter.count(m) || full.count(m);
  }

  bool only_full() const
  {
    return early_ender.empty();
  }

  Match * early()
  {
    return !early_ender.empty() ? *early_ender.begin() : *full.begin();
  }
  
  Match * late()
  {
    return !late_starter.empty() ? *late_starter.begin() : *full.begin();
  }

  // We suppose that m \in this
  bool is_early_or_late(Match* m)
  {
    return only_full()
      || early()->projected_end(t)==m->projected_end(t)
      || late()->projected_start(t)==m->projected_start(t);
  }

  
  // for debug
  void display() const
  {
    cout << "Early\n";
    for(auto m : early_ender)
      cout << *m;
    cout << "Late\n";
    for(auto m : late_starter)
      cout << *m;
    cout << "Full\n";
    for(auto m : full)
      cout <<*m;

  }
  // for debug
  void display_adress() const
  {
    cout << "Early " << early_ender.size() <<"\n";
    for(auto m : early_ender)
      cout << m << " ";
    cout << "\nLate " << late_starter.size()<< "\n";
    for(auto m : late_starter)
      cout << m << " ";
    cout << "Full\n";
    for(auto m : full)
      cout << m << " ";
  }
  
};






