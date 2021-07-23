#pragma once

#include <iostream>
#include <set>
#include <vector>


using namespace std;

class Contig
{
  friend ostream& operator<<(ostream&, const Contig&);
private:
  vector<unsigned> str;
  const size_t contig_size;
  const string name;
		
public:

  // indicates if this contig belongs to the set T
  const bool t;

  
  Contig(bool t, string &&name,size_t size=0) : t(t), name(move(name)), contig_size(size) {}

  Contig(vector<unsigned>&& arr, string &&name, bool t) : str(move(arr)), contig_size(arr.size()), name(move(name)), t(t) {}

  bool operator<(const Contig &c) const
  {
    return name<c.name;
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

  /*int& atrev(int i)
    {
    return str[str.size() - 1 - i];	//backcheck this
    }*/


  size_t size() const
  {
    return contig_size;
  }



};
typedef set<Contig> AssemblySet;


ostream &operator<<(ostream& os, const Contig& contig)
{
  for(int c : contig.str)
    os<<(char)c;
  return os;
}
