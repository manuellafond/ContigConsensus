#pragma once

#include <array>
#include <bitset>
#include <iostream>

using namespace std;


  class Nuc : public array<unsigned,4>
  {
    friend Nuc operator+(const Nuc&, const Nuc &);
    
    Nuc() : array<unsigned int, 4>{{0,0,0,0}} {}
    
    public:
    Nuc(char c) : Nuc()
    {
	*this=c; 
    }

    void operator=(char c)
      {
	switch (c) {
	case 'N':
	  for(unsigned i=0; i<4;++i)
	    (*this)[i]++;
	  break;
	case 'A':
	  (*this)[0]++;
	  break;
	case 'T':
	  (*this)[1]++;
	  break;
	case 'C':
	  (*this)[2]++;
	  break;
	case 'G':
	  (*this)[3]++;
	  break;
	}
      }

    bool operator==(const Nuc &n)
      {
	for(size_t i =0; i<4;++i)
	  if((*this)[i] && n[i])
	    return true;
	return false;
      }

    void operator+=(const Nuc&n)
    {
      for(size_t i =0; i<4;++i)
	(*this)[i]+=n[i];
    }

    operator char() const
    {
      unsigned m=0;
      for(unsigned i =1; i <4;i++)
	if((*this)[i]>(*this)[m])
	  m=i;
      switch (m) {
      case 0:
	return 'A';
	break;
      case 1:
	return 'T';
	break;
      case 2:
	return 'C';
	break;
      default:
	return 'G';
	break;
      }
      
    }

   
  };               

   Nuc operator+(const Nuc&n1, const Nuc &n2)
    {
      Nuc n;
      for(size_t i =0; i<4;++i)
	n[i]=n1[i]+n2[i];
      return n;
    }
