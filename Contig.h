#pragma once

#include <iostream>
#include <vector>


using namespace std;

class Contig
{
  friend ostream& operator<<(ostream&, const Contig&);
private:
  vector<unsigned> str;
		
public:

  // indicates to which set of contigs its belongs
  const unsigned t; 
  

  Contig(unsigned t, size_t size=0) : t(t)
  {
    str.resize(size);
  }

  Contig(vector<unsigned>&& arr, unsigned t) : t(t)
  {
    str = move(arr);
  }


  void resize(size_t dim) 
  {
    str.resize(dim);
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


  size_t size()
  {
    return str.size();
  }



};
typedef vector<Contig> AssemblySet;


ostream &operator<<(ostream& os, const Contig& contig)
{
  for(int c : contig.str)
    os<<(char)c;
  return os;
}
