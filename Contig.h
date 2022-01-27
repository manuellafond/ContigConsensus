#pragma once

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <tuple>
#include <vector>

#include "Nucleotide.h"

using namespace std;

class Contig
{
  friend ostream& operator<<(ostream&, const Contig&);
private:

  vector<Nuc> sequence;
  const string name;
  const unsigned set_id;
  size_t sequence_size;


  struct Component{
    string name;
    unsigned set_id;
    bool is_reversed;
    size_t shift, contig_size;

    Component(size_t contig_size,string name, unsigned set_id, size_t shift,bool is_reversed=false) :  name(name), set_id(set_id), is_reversed(is_reversed), shift(shift), contig_size(contig_size) {}

    //for debug
    void display()
    {
      cout << "{" << name << "," << set_id << "," << is_reversed << "," << shift << "}" << endl;
    }
  };
  vector<Component> component_contigs;

  
public:
  

  Contig(unsigned set_id, string &&name,vector<Nuc> && sequence) : name(name), set_id(set_id), sequence(move(sequence)) {
    sequence_size=this->sequence.size();
    component_contigs.emplace_back(sequence_size,name,set_id,0);
  }

    // merge two contigs
  Contig(Contig && c1, Contig && c2, bool is_c1_reversed, bool is_c2_reversed, size_t overlapping_size, unsigned set_id) :  name(c1.name + "|" + c2.name), set_id(set_id), sequence_size(c1.size())
  {
    if(!is_c1_reversed){
      sequence=move(c1.sequence);
      component_contigs = move(c1.component_contigs);
    }
    else {
      for(unsigned i=c1.size();i>0;--i)
	sequence.push_back(c1.sequence[i-1]);
           for(auto &C : c1.component_contigs)
	     component_contigs.emplace_back(C.contig_size,C.name,C.set_id,c1.sequence_size-1-C.shift,!C.is_reversed);
    }


    unsigned shift=c1.size()-overlapping_size;
    if(!is_c2_reversed)
       for(auto &C : c2.component_contigs)
	 component_contigs.emplace_back(C.contig_size,C.name,C.set_id,shift+C.shift,C.is_reversed);
    else
      for(auto &C : c2.component_contigs)
	component_contigs.emplace_back(C.contig_size,C.name,C.set_id,shift+c2.sequence_size-1-C.shift,!C.is_reversed);

    if(is_c2_reversed){
      for(size_t i=0;i<overlapping_size;++i){
	if(c2.size()<=i) return;
	sequence[sequence.size()-overlapping_size+i]+=c2.sequence[c2.sequence.size()-1-i];
	// cout << sequence.size()-overlapping_size+i << "|"  << c2.sequence.size()-1-i  <<":";
	// display_sequence();
      }
      for(size_t i=overlapping_size;i<c2.sequence.size();++i)
	sequence.push_back(move(c2.sequence[c2.sequence.size()-1-i]));

    }
    else { 
      for(size_t i=0;i<overlapping_size;++i){
	if(c2.size()<=i) return;
	sequence[sequence.size()-overlapping_size+i]+=c2.sequence[i];
	// cout << i << "|"  << sequence.size()-overlapping_size+i <<":";
	// display_sequence();
      }
      
      for(size_t i=overlapping_size;i<c2.sequence.size();++i)
	sequence.push_back(move(c2.sequence[i]));
      
    }
    sequence_size=sequence.size();
  }



  unsigned get_set_id() const
  {
    return set_id;
  }
  
  bool operator<(const Contig &c) const
  {
    return name<c.getName();
  }
  

  const string & getName() const
  { 
    return name;
  } 


  char at(size_t i, bool isReverse = false)
  {
    if (!isReverse)
      return sequence[i];
    else
      return sequence[sequence.size() - 1 - i];
  }

  char getNuc(size_t i)
  {
    return sequence[i];
  }


  size_t size() const
  {
    return sequence_size;
  }

  void display_sequence() const // for debug
  {
    for(auto &N : this->sequence)
      cout << (char)N;
    cout << endl;
  }
  void display_component() const // for debug
  {
    for(auto c : this->component_contigs)
      c.display();
  }

  vector<Component>& getComponent() {
    return component_contigs;
  }
  

};

bool operator<(const unique_ptr<Contig> &c, const string &s)
{
  return c->getName() < s;
}
bool operator<(const string &s,const unique_ptr<Contig> &c)
{
  return s< c->getName();
}

typedef map<unsigned, set<unique_ptr<Contig>,less<>>> AssemblySet;




inline ostream &operator<<(ostream& os, const Contig& contig)
{
  os << ">"<<contig.name << endl;
  for(Nuc c : contig.sequence)
    os<<(char)c;
  os << endl;
  return os;
}
