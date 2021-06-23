#pragma once

#include <dirent.h>
#include <iostream>
#include <fstream>
#include <limits>

#include "Contig.h"

#define DEBUG true

using namespace std;




vector<Contig> parseFile(const char* fileName, int t)
{
	ifstream file(fileName);
	if(!file)
	{
		cout << "Error: could not open " << fileName << endl;
		exit(EXIT_FAILURE);
	}
  
	vector<Contig> contigs;
	while(!file.eof())
	{
		vector<unsigned> seq;
       
		// ignore commentary line
		file.ignore(numeric_limits<streamsize>::max(),'\n');
		
		int c;
		do {
			c = file.get();
			if(c == '\n' || c == EOF|| c == '>' || c == '*') continue;
			seq.push_back(c);
			
		} while (c!='>' && c != EOF);
		if(seq.size()>0)
		    contigs.push_back(Contig(move(seq),t));
		
	}

	#if DEBUG==true
	for(const Contig &cont : contigs)
	  cout << cont << "\n";
	#endif

	return contigs;
}

// vector<AssemblySet> treatDirectory(const char *dirName)
// {
// 	vector<AssemblySet> assemblySet;
// 	DIR* dir;
// 	struct dirent *ent;
// 	dir=opendir(dirName);
//         if (dir == NULL) {
// 		cout << "Error: could not open directory " << dirName << "\n";
// 		exit(EXIT_FAILURE);
//         }

// 	while((ent=readdir(dir))!=NULL){
// 		if(ent->d_type!=DT_REG)
// 			continue;
		
// 		string f_name = ent->d_name;
// 		if(f_name.find(".fasta")== string::npos
// 		   || f_name.substr(f_name.find_last_of(".")) != ".fasta")
// 			continue;
		
// 		cout << "Read fasta file: " << f_name <<"\n";
// 		string file = dirName;
// 		file.append("/");
// 		file.append(ent->d_name);
// 		assemblySet.push_back(parseFile(file.c_str()));
// 	}
// 	return assemblySet;
// }
