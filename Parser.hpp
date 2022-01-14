#pragma once

#include <dirent.h>
#include <iostream>
#include <fstream>
#include <limits>
#include <memory>
#include <set>
#include <stdexcept>
#include <unordered_set>
#include <utility>
#include <vector>

#include "Contig.h"
#include "Match.h"

#define DEBUG false

using namespace std;

vector<Contig> parseFile(const char* fileName, unsigned id_set, AssemblySet &assembly_sets)
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


void parseMatches(const char * fileName, AssemblySet & assembly_sets, MatchMatrix & matches, unsigned id_set1=1, unsigned id_set2=2)
{
  ifstream file(fileName);
  if(!file) {
      cout << "Error: could not open " << fileName << endl;
      exit(EXIT_FAILURE);
  }

  // insert keys
  auto &v = matches[id_set1][id_set2];
  
  size_t score, startS, endS, lengthS, startT, endT, lengthT;
  string nameS, nameT;
  while (file >> score) {

    file >> nameS;
    file >> startS;
    file >> endS;
    file >> lengthS;
    file >> nameT;
    file >> startT;
    file >> endT;
    file >> lengthT;

    Contig *s, *t;
    auto a = assembly_sets[id_set1].find(nameS);
    if(a!=assembly_sets[id_set1].end())
      s=a->get();
    else s=assembly_sets[id_set1].insert(move(make_unique<Contig>(id_set1,move(nameS),lengthS))).first->get();

    a = assembly_sets[id_set2].find(nameT);
    if(a!=assembly_sets[id_set2].end())
      t=a->get();
    else t=assembly_sets[id_set2].emplace(make_unique<Contig>(id_set2,move(nameT),lengthT)).first->get();

    
    // // modify the match length so that it is a suffix/prefix match
    // size_t resize_l = min({startS,endS,startT,endT})-1,
    //   resize_r = min({lengthS-startS,lengthS-endS,lengthT-startT,lengthT-endT});
    
    // if(startS<endS){
    //   startS-=resize_l;
    //   endS+=resize_r;
    // } else{
    //   startS+=resize_r;
    //   endS-=resize_l;
    // }

    // if(startT<endT){
    //   startT-=resize_l;
    //   endT+=resize_r;
    // } else{
    //   startT+=resize_r;
    //   endT-=resize_l;
    // }

    v.emplace_back(s,t,startS,endS,startT,endT,score);
    
    
  }

  cout << "test\n";
  for(auto &m : v)
    m.display_contig_names();
  
  cout << "Size: " << assembly_sets[id_set1].size() << " " << assembly_sets[id_set2].size() << " " << v.size() <<  endl;

}

void output(const char * fileName, vector<const Match*>& matches)
{

  ofstream file(fileName);
  if(!file) {
    cout << "Error: could not open " << fileName << endl;
    exit(EXIT_FAILURE);
  }
  //  file << "#" << fileTName << " vs " << fileSName << endl;

  for(size_t i = 0; i<matches.size(); ++i){
    file << "#Match" << i+1<<"-score-" << matches[i]->score<<endl;
    for(unsigned t = 0; t<=1; ++t){
      //file << "#" << (t==0 ? fileTName : fileSName) << endl;
      file <<">" << matches[i]->contig(t)->getName() << ": " 
	   << matches[i]->start(t) << " " << matches[i]->end(t) 
	    <<endl;

    }
    file << endl;
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
