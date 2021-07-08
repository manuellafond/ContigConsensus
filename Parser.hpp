#pragma once

#include <dirent.h>
#include <iostream>
#include <fstream>
#include <limits>

#include "Match.h"

#define DEBUG false

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
	
	file.get(); // do not read the first '>'
	while(!file.eof())
	{
		vector<unsigned> seq;
       
		// ignore commentary line
		//		file.ignore(numeric_limits<streamsize>::max(),'\n');
		string name;
		getline(file,name);
		
		int c;
		do {
			c = file.get();
			if(c == '\n' || c == EOF|| c == '>' || c == '*') continue;
			if((char)c != 'A' && (char)c != 'T' && (char)c != (char)'C' && c != (char)'G') continue;
			seq.push_back(c);
			
		} while (c!='>' && c != EOF);
		if(seq.size()>0)
		  contigs.push_back(Contig(move(seq),move(name),t));
		
	}

	#if DEBUG==true
	for(const Contig &cont : contigs)
	  cout << cont << "\n";
	#endif

	return contigs;
}

void output(const char * fileName, const char * fileTName, const char * fileSName, vector<Match*>& matches)
{

  ofstream file(fileName);
  if(!file) {
    cout << "Error: could not open " << fileName << endl;
    exit(EXIT_FAILURE);
  }
  file << "#" << fileTName << " vs " << fileSName << endl;

  for(size_t i = 0; i<matches.size(); ++i){
    file << "#Match" << i+1<<"-score-" << matches[i]->score<<endl;
    for(unsigned t = 0; t<=1; ++t){
      file << "#" << (t==0 ? fileTName : fileSName) << endl;
      file <<">" << matches[i]->contig(t)->getName() << ": " 
	   << matches[i]->start(t) << " " << matches[i]->end(t) ;
      
      for (size_t x = 0; x < matches[i]->length; ++x) {
	if (x % 80 == 0) 
	  file << endl;
	if((char)matches[i]->contig(t)->at(matches[i]->start(t) + x, matches[i]->is_reverse(t))=='\n')
	  cout << "c";
	file << (char)matches[i]->contig(t)->at(matches[i]->start(t) + x, matches[i]->is_reverse(t));
      }
      file << endl <<endl;

    }
  }
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
