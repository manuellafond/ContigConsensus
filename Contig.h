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
  
	Contig() : Contig(0)
	{
		
	}

	Contig(int size)
	{
		str.resize(size);
	}

	Contig(vector<unsigned>&& arr)
	{
		str = move(arr);
	}


	void resize(int dim)
	{
		str.resize(dim);
	}


	unsigned& at(int i, bool isReverse = false)
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


	int size()
	{
		return str.size();
	}



};
typedef vector<Contig> AssemblySet;
