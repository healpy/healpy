/*
 *  This file is part of libcxxsupport.
 *
 *  libcxxsupport is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  libcxxsupport is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with libcxxsupport; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/*
 *  libcxxsupport is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*
 *  This file contains functionality related to wall-clock timers
 *
 *  Copyright (C) 2010, 2011, 2012 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include <iostream>
#include <utility>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include "walltimer.h"
#include "walltime_c.h"
#include "error_handling.h"

using namespace std;

void wallTimer::start()
  { start(wallTime()); }
void wallTimer::stop()
  { stop(wallTime()); }
double wallTimer::acc() const
  { return acc(wallTime()); }

int wallTimerSet::getIndex(const string &name)
  {
  maptype::const_iterator it = lut.find(name);
  if (it!=lut.end())
    return it->second;
  timer.push_back(wallTimer());
  lut[name]=timer.size()-1;
  return timer.size()-1;
  }

void wallTimerSet::start(int index)
  { timer[index].start(); }
void wallTimerSet::stop(int index)
  { timer[index].stop(); }
void wallTimerSet::stopstart(int index1, int index2)
  { double t=wallTime(); timer[index1].stop(t); timer[index2].start(t); }
void wallTimerSet::reset(int index)
  { timer[index].reset(); }
double wallTimerSet::acc(int index)
  { return timer[index].acc(); }
void wallTimerSet::start(const string &name)
  { start(getIndex(name)); }
void wallTimerSet::stop(const string &name)
  { stop(getIndex(name)); }
void wallTimerSet::stopstart(const string &name1, const string &name2)
  { stopstart(getIndex(name1),getIndex(name2)); }
void wallTimerSet::reset(const string &name)
  { reset(getIndex(name)); }
double wallTimerSet::acc(const string &name)
  { return acc(getIndex(name)); }

void wallTimerSet::report() const
  {
  cout << "\nWall clock timer report:" << endl;
  for (maptype::const_iterator it=lut.begin(); it!=lut.end(); ++it)
    printf("  %-15s: %10.5fs\n", it->first.c_str(), timer[it->second].acc());
  cout << "End wall clock timer report\n" << endl;
  }

wallTimerSet wallTimers;

namespace {

class tstack_node;

typedef map<string,tstack_node>::iterator Ti;
typedef map<string,tstack_node>::const_iterator Tci;
typedef pair<Tci,double> Tipair;

class tstack_node
  {
  public:
    tstack_node *parent;
    wallTimer wt;
    string name;
    map<string,tstack_node> child;

    tstack_node(const string &name_, tstack_node *parent_)
      : parent(parent_), name(name_) {}

    int max_namelen() const
      {
      int res=name.length();
      for (Tci it=child.begin(); it!=child.end(); ++it)
        res=max(res,it->second.max_namelen());
      return res;
      }
  };

tstack_node tstack_root("root",0);
tstack_node *curnode=0;
double overhead=0.;

struct timecomp
  {
  bool operator() (const Tipair &a, const Tipair &b) const
    { return a.second>b.second; }
  };

void tstack_report(const tstack_node &node, const string &indent, int twidth,
  int slen)
  {
  double total=node.wt.acc();
  vector<Tipair> tmp;
  for (Tci it=node.child.begin(); it!=node.child.end(); ++it)
    tmp.push_back(make_pair(it,it->second.wt.acc()));

  if (tmp.size()>0)
    {
    sort(tmp.begin(),tmp.end(),timecomp());
    double tsum=0;
    printf("%s|\n", indent.c_str());
    for (unsigned i=0; i<tmp.size(); ++i)
      {
      printf("%s+- %-*s:%6.2f%% (%*.4fs)\n",indent.c_str(),slen,
        (tmp[i].first->first).c_str(), 100*tmp[i].second/total,twidth,
        tmp[i].second);
      tstack_report(tmp[i].first->second,indent+"|  ",twidth,slen);
      tsum+=tmp[i].second;
      }
    printf("%s+- %-*s:%6.2f%% (%*.4fs)\n%s\n",indent.c_str(),slen,
      "<unaccounted>",100*(total-tsum)/total,twidth,total-tsum,indent.c_str());
    }
  }

} // unnamed namespace

void tstack_push(const string &name)
  {
  double t0=wallTime();
  if (curnode==0) curnode=&tstack_root;
  Ti it=curnode->child.find(name);
  if (it==curnode->child.end())
    it=curnode->child.insert (make_pair(name,tstack_node(name,curnode))).first;
  curnode=&(it->second);
  double t1=wallTime();
  curnode->wt.start(0.5*(t0+t1));
  overhead+=t1-t0;
  }
void tstack_pop(const string &name)
  {
  double t0=wallTime();
  planck_assert(curnode && (curnode->name==name), "invalid tstack operation");
  double t1=wallTime();
  curnode->wt.stop(0.5*(t0+t1));
  curnode=curnode->parent;
  overhead+=t1-t0;
  }
void tstack_pop()
  {
  double t0=wallTime();
  planck_assert(curnode, "invalid tstack operation");
  double t1=wallTime();
  curnode->wt.stop(0.5*(t0+t1));
  curnode=curnode->parent;
  overhead+=t1-t0;
  }
void tstack_replace(const string &name2)
  {
  double t0=wallTime();
  planck_assert(curnode, "invalid tstack operation");
  tstack_node *savenode=curnode;
  curnode=curnode->parent;
  Ti it=curnode->child.find(name2);
  if (it==curnode->child.end())
    it=curnode->child.insert(make_pair(name2,tstack_node(name2,curnode))).first;
  curnode=&(it->second);
  double t1=wallTime();
  double t=0.5*(t0+t1);
  savenode->wt.stop(t);
  curnode->wt.start(t);
  overhead+=t1-t0;
  }
void tstack_replace(const string &name1, const string &name2)
  {
  planck_assert(curnode && (curnode->name==name1), "invalid tstack operation");
  tstack_replace(name2);
  }

void tstack_report(const string &stem)
  {
  const tstack_node *ptr = 0;
  for (Tci it=tstack_root.child.begin(); it!=tstack_root.child.end(); ++it)
    if (it->first==stem) ptr=&(it->second);
  planck_assert(ptr,"invalid stem");
  int slen=string("<unaccounted>").size();
  slen = max(slen,ptr->max_namelen());

  double total=ptr->wt.acc();
  printf("\nTotal wall clock time for '%s': %1.4fs\n",stem.c_str(),total);

  int logtime=max(1,int(log10(total)+1));
  tstack_report(*ptr,"",logtime+5,slen);

  printf("\nAccumulated timing overhead: approx. %1.4fs\n",overhead);
  }
