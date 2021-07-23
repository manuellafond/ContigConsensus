#pragma once

#include <algorithm>
#include<bits/stdc++.h> 
#include <set>

#include "Match.h"

#define DEBUG_HUNG false

struct data {
  unsigned nb_selected_zero;
  vector<vector<pair<int, const Match*>>> mat;
  vector<int> selected_zero_r, selected_zero_c, zero_prime_r, zero_prime_c;
  vector<bool> mark_r, mark_c;

  data(int max_value,size_t size ) : nb_selected_zero(0),
				     mat(size, vector<pair<int,const Match*>>(size,make_pair(max_value, nullptr))),
				     selected_zero_r(size,-1), selected_zero_c(size,-1),
				     zero_prime_r(size,-1), zero_prime_c(size,-1),
				     mark_r(size,false), mark_c(size,false) {}

  void display()
  {
    // for debug
    cout << "\nNb: " <<nb_selected_zero <<endl;
    cout << "  ";
    for(bool m : mark_c)
      if(m)
	cout << " M ";
      else cout << "   ";
    cout <<endl;
    cout << "  ";
    for(auto m : mark_c)
      cout << "+++";
    cout << endl;
    for(int i =0; i <mat.size(); ++i){
      if(mark_r[i])
	cout << "M";
      else cout << " ";
      cout << "+";
      for(int j =0; j <mat.size();++j){
	if(selected_zero_r[i]==j)
	  cout <<"|";
	else cout << " ";
	cout << mat[i][j].first;
	if(selected_zero_r[i]==j)
	  cout <<"|";
	else if(zero_prime_r[i]==j)
	  cout <<"'";
	else cout << " ";
      }
      cout << endl;
    }
  }
 
};





bool select_zero(data &d)
{
  //d.display();
  unsigned nb_zero=0;
  for(int i=0; i<d.mat.size();++i){
    if(d.selected_zero_r[i]!=-1)
      continue;
    for(int j=0; j<d.mat.size();++j)
      if(d.mat[i][j].first==0 && d.selected_zero_c[j]==-1){
	d.selected_zero_r[i]=j;
	d.selected_zero_c[j]=i;
	nb_zero++;
	break;
      }
    }
  d.nb_selected_zero+=nb_zero;

  //d.display();

  if(d.nb_selected_zero==d.mat.size())
    return true; // optimality!
  if(nb_zero>0){
    for(auto *v : {&d.zero_prime_c,&d.zero_prime_r})
      v->assign(d.mat.size(),-1);
    for(auto *v : {&d.mark_c,&d.mark_r})
      v->assign(d.mat.size(),false);
  }

  
  return false;
}


void augment(data &d,int i, int j)
{
  //  cout << "augment\n";
  // step 2'
  vector<array<int,2>> invert_s,invert_p;
  invert_p.push_back({i,j});
  while(d.selected_zero_c[j]!=-1){
    i=d.selected_zero_c[j];
    invert_s.push_back({i,j});
    j=d.zero_prime_r[i];
    invert_p.push_back({i,j});
  }
  // not sure the first loop is usefull
  for(auto z : invert_s){
    d.selected_zero_r[z[0]]=-1;
    d.selected_zero_c[z[1]]=-1;
  }
  for(auto z : invert_p){
    d.selected_zero_r[z[0]]=z[1];
    d.selected_zero_c[z[1]]=z[0];
  }

  d.nb_selected_zero++;
}

bool prime(data &d)
{
  for(int i=0;i<d.mat.size();++i){
    if(d.selected_zero_c[i]!=-1)
      d.mark_c[i]=true;
    if(d.selected_zero_r[i]!=-1)
      d.mark_r[i]=false;
  }

  bool unmarked_zero = false;
  while(!unmarked_zero){
    unmarked_zero=true;
  for(int i=0; i<d.mat.size();++i)
    for(int j =0; j<d.mat.size();++j){
      if (d.mat[i][j].first == 0 && !d.mark_r[i] && !d.mark_c[j]) {
	unmarked_zero=false;
	
        d.zero_prime_r[i]=j;
	d.zero_prime_c[j]=i;

        if (d.selected_zero_r[i] == -1) {
          augment(d, i, j);
          return false;
        }
        d.mark_c[d.selected_zero_r[i]] = false;

        d.mark_r[i] = true;
      }
    }
  }
  return true;
}

void substract_min(data &d)
{
  int min_tmp=INT_MAX;
  for (int i = 0; i < d.mat.size(); ++i) {
    if(d.mark_r[i]) continue;
    for (int j = 0; j < d.mat.size(); ++j) {
      if(d.mark_c[j]) continue;
      min_tmp=min(min_tmp,d.mat[i][j].first);
    }
  }



  for (int i = 0; i < d.mat.size(); ++i) {
    for (int j = 0; j < d.mat.size(); ++j) {
      if(d.mark_r[i])
	d.mat[i][j].first+=min_tmp;
      if(!d.mark_c[j])
	d.mat[i][j].first-=min_tmp;
    }
  }

}


std::vector<const Match*> hungarian_algorithm(AssemblySet &T, AssemblySet &S, const Match::MM_map& matches, int max_value, unsigned &score)
{

  #if DEBUG_HUNG
  cout << "S: ";
  for(auto &s : S)
    cout << s << " ";
  cout <<endl;
  cout << "T: ";
  for(auto &t : T)
    cout << t << " ";
  cout <<endl;
  #endif


  data d(max_value,max(S.size(),T.size()));

  
  size_t i=0;
  for(const Contig &s : S){
    size_t j=0;
    for(const Contig &t : T) {
      try {
	d.mat[i][j].first -= matches.at(&s).at(&t)->score;
	d.mat[i][j].second = matches.at(&s).at(&t);

      } catch (out_of_range) {
	//        d.mat[i][j].first = 0;
      }
      ++j;
    }
    ++i;
  }
  for (int i = 0; i < d.mat.size(); ++i) {
    int min_tmp= INT_MAX;
    for (int j = 0; j < d.mat.size(); ++j)
      min_tmp=min(min_tmp,d.mat[i][j].first);
    for (int j = 0; j < d.mat.size(); ++j)
      d.mat[i][j].first-=min_tmp;
  }

  for (int j = 0; j < d.mat.size(); ++j) {
    int min_tmp=INT_MAX;
    for (int i = 0; i < d.mat.size(); ++i)
      min_tmp=min(min_tmp,d.mat[i][j].first);
    for (int i = 0; i < d.mat.size(); ++i)
      d.mat[i][j].first-=min_tmp;
  }
  while (!select_zero(d)) {
    if (prime(d))
      substract_min(d);
  }

#if DEBUG_HUNG
  cout << "Selected matches:\n";
  for (int i = 0; i < d.mat.size(); ++i) {
    if(i>S.size()-1)
      cout << "(empty,";
    else cout <<"(" << S[i] <<",";
    if(d.selected_zero_r[i]>T.size()-1)
       cout << "empty)";
    else cout << T[d.selected_zero_r[i]] << ")";
    cout << ": score -> ";
    if(i>S.size()-1 || d.selected_zero_r[i]>T.size()-1)
      cout << "0\n";
    else
      try {
        cout << matches.at(&S[i]).at(&T[d.selected_zero_r[i]])->score << endl;
      }
      catch (out_of_range) {
	cout << "0\n";
      }
  }
#endif

  
  vector<const Match*> selected_matches;
  for (size_t i = 0; i < d.mat.size(); ++i)
        if (d.mat[i][d.selected_zero_r[i]].second!= nullptr) {
	  selected_matches.push_back(d.mat[i][d.selected_zero_r[i]].second);
	  score+=d.mat[i][d.selected_zero_r[i]].second->score;
        }
      
    
  cout << "size:" << selected_matches.size() << endl;
  return selected_matches;
}




