#pragma once
#include <iostream>
#include <fstream>
#include <limits>
#include <map>
#include <string.h>

using namespace std;

class CostMap
{
	map<unsigned, map<unsigned, int>> cost;
public:
	CostMap() {}
	CostMap(const char* fileName){
		if(strcmp(fileName, "-1")==0)
			return;
		ifstream file(fileName);
		if(!file)
		{
			cout << "Error: could not open " << fileName << endl;
			exit(EXIT_FAILURE);
		}
		while(!file.eof()){
			char c1,c2;
			int value;
			file>>c1>>c2>>value;
			file.ignore(numeric_limits<streamsize>::max(),'\n');
			cost[c1][c2]=value;
		}
	}

	int operator() (unsigned c1, unsigned c2){
		if(cost.find(c1)==cost.end() || cost.find(c2)==cost.end())
			return c1==c2;
		return cost[c1][c2];
	}
       
};
