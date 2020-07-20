#pragma once

#include <vector>

using namespace std;

class Contig
{
private:
	vector<int> str;

public:

	Contig() : Contig(0)
	{
		
	}

	Contig(int size)
	{
		str.resize(size);
	}

	Contig(vector<int> arr)
	{
		str = arr;
	}


	void resize(int dim)
	{
		str.resize(dim);
	}


	int& at(int i, bool isReverse = false)
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

