#pragma once

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <tuple>
#include <vector>


using namespace std;

class Contig
{
  friend ostream& operator<<(ostream&, const Contig&);
private:

  const size_t contig_size;
  const string name;
  const unsigned set_id;
  
  struct Component{
    size_t contig_size;
    string name;
    unsigned set_id;
    bool is_reversed;
    size_t shift;

    Component(size_t contig_size, string name, unsigned set_id, size_t shift,bool is_reversed=false) : contig_size(contig_size), name(name), set_id(set_id), is_reversed(is_reversed), shift(shift) {}

    //for debug
    void display()
    {
      cout << "{" << contig_size << "," << name << "," << set_id << "," << is_reversed << "," << shift << "}" << endl;
    }
  };
  vector<Component> component_contigs;
  
  
  vector<unsigned> str;

public:

  
  
  Contig(unsigned set_id, string &&name,size_t size=0) : contig_size(size), name(name), set_id(set_id), component_contigs{{Component(size,name,set_id,0)}}  {}

  // merge two contigs
  Contig(Contig && c1, Contig && c2, bool is_c2_reversed, size_t overlapping_size, string name, unsigned set_id) : contig_size(c1.size()+c2.size()-overlapping_size),name(name), set_id(set_id), component_contigs(move(c1.component_contigs))
  {
    size_t shift=c1.size()-overlapping_size;
    if(is_c2_reversed){
      for(size_t i=c2.component_contigs.size();i>0;--i){
	auto & cc = c2.component_contigs[i-1];
	cc.shift=shift+c2.size()-(cc.shift+cc.contig_size);
	cc.is_reversed^=true;
	component_contigs.push_back(move(cc));
      }

    }
    else {
      for(size_t i=0;i<c2.component_contigs.size();++i){
	c2.component_contigs[i].shift+=shift;
	component_contigs.push_back(move(c2.component_contigs[i]));
      }
    }
    sort(component_contigs.begin(),component_contigs.end(),[](const Component &c1, const Component &c2){return c1.shift<c2.shift;});

    display_component();
  }

  unsigned get_set_id() const
  {
    return set_id;
  }
  
  bool operator<(const Contig &c) const
  {
    return name<c.getName();
  }
  
  void resize(size_t dim) 
  {
    str.resize(dim);
  }

  const string & getName() const
  { 
    return name;
  } 


  unsigned& at(size_t i, bool isReverse = false)
  {
    if (!isReverse)
      return str[i];
    else
      return str[str.size() - 1 - i];
  }


  size_t size() const
  {
    return contig_size;
  }

  //for debug
  void display_component() const
  {
    for(auto c : this->component_contigs)
      c.display();
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
  for(int c : contig.str)
    os<<(char)c;
  return os;
}
