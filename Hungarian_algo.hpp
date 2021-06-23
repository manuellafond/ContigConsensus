#pragma once

#include <algorithm>
#include<bits/stdc++.h> 
#include <set>

#include "Match.h"

#define DEBUG_HUNG true

struct data {
  unsigned nb_selected_zero;
  vector<vector<int>> mat;
  vector<int> selected_zero_r, selected_zero_c, zero_prime_r, zero_prime_c;
  vector<bool> mark_r, mark_c;

  data(int max_value,size_t size ) : nb_selected_zero(0),
				     mat(size, vector<int>(size,max_value)),
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
	cout << mat[i][j];
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

  unsigned nb_zero=0;
  for(int i=0; i<d.mat.size();++i){
    if(d.selected_zero_r[i]!=-1)
      continue;
    for(int j=0; j<d.mat.size();++j)
      if(d.mat[i][j]==0 && d.selected_zero_c[j]==-1){
	d.selected_zero_r[i]=j;
	d.selected_zero_c[j]=i;
	nb_zero++;
	break;
      }
    }
  d.nb_selected_zero+=nb_zero;



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

  
  for(int i=0; i<d.mat.size();++i)
    for(int j =0; j<d.mat.size();++j){
      if (d.mat[i][j] == 0 && d.selected_zero_r[i] != j) {
        d.zero_prime_r[i]=j;
	d.zero_prime_c[j]=i;

        if (d.selected_zero_r[i] == -1) {
	  augment(d,i,j);
	  return false;
        }
	d.mark_c[d.selected_zero_r[i]]=false;
	d.mark_r[i]=true;
	
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
      min_tmp=min(min_tmp,d.mat[i][j]);
    }
  }

  for (int i = 0; i < d.mat.size(); ++i) {
    for (int j = 0; j < d.mat.size(); ++j) {
      if(d.mark_r[i])
	d.mat[i][j]+=min_tmp;
      if(!d.mark_c[j])
	d.mat[i][j]-=min_tmp;
    }
  }
  
}


std::vector<Match*> hungarian_algorithm(AssemblySet &T, AssemblySet &S, MM_map& matches, int max_value)
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

  for(int i = 0; i<S.size();++i)
    for (int j = 0; j < T.size(); ++j)
      if(matches[&S[i]][&T[j]]!=nullptr)
	d.mat[i][j] -= matches[&S[i]][&T[j]]->score;
    

  for (int i = 0; i < d.mat.size(); ++i) {
    int min_tmp=INT_MAX;
    for (int j = 0; j < d.mat.size(); ++j)
      min_tmp=min(min_tmp,d.mat[i][j]);
    for (int j = 0; j < d.mat.size(); ++j)
      d.mat[i][j]-=min_tmp;
  }

  for (int j = 0; j < d.mat.size(); ++j) {
    int min_tmp=INT_MAX;
    for (int i = 0; i < d.mat.size(); ++i)
      min_tmp=min(min_tmp,d.mat[i][j]);
    for (int i = 0; i < d.mat.size(); ++i)
      d.mat[i][j]-=min_tmp;
  }
  
  while(!select_zero(d))
    if(prime(d))
      substract_min(d);

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
    if(i>S.size()-1 || d.selected_zero_r[i]>T.size()-1 || matches[&S[i]][&T[d.selected_zero_r[i]]]==nullptr)
      cout << "0\n";
    else cout << matches[&S[i]][&T[d.selected_zero_r[i]]]->score<< endl;

  }
#endif
  
  vector<Match*> selected_matches;
  for (int i = 0; i < d.mat.size(); ++i)
    if (i < S.size() && d.selected_zero_r[i] < T.size()) {
      if(matches[&S[i]][&T[d.selected_zero_r[i]]] != nullptr)
	selected_matches.push_back(matches[&S[i]][&T[d.selected_zero_r[i]]]);
    } 

  return selected_matches;
}
