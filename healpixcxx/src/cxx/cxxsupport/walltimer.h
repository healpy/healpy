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

/*! \file walltimer.h
 *  Functionality related to wall-clock timers
 *
 *  Copyright (C) 2010, 2011 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef PLANCK_WALLTIMER_H
#define PLANCK_WALLTIMER_H

#include <string>
#include <map>
#include <vector>

class wallTimer
  {
  private:
    double t_acc, t_started;
    bool running;

  public:
    wallTimer() : t_acc(0.), t_started(0.), running(false) {}
    void start(double wtime_now)
      { if (!running) { t_started=wtime_now; running=true; } }
    void start();
    void stop(double wtime_now)
      { if (running) { t_acc+=wtime_now-t_started; running=false; } }
    void stop();
    void reset() { t_acc=t_started=0.; running=false;}
    double acc(double wtime_now) const
      { return running ? t_acc+wtime_now-t_started : t_acc; }
    double acc() const;
  };

class wallTimerSet
  {
  public:
    typedef std::map<std::string,int> maptype;

  private:
    maptype lut;
    std::vector<wallTimer> timer;

  public:
    int getIndex(const std::string &name);
    void start(int index);
    void stop(int index);
    void stopstart(int index1, int index2);
    void reset(int index);
    double acc(int index);
    void start(const std::string &name);
    void stop(const std::string &name);
    void stopstart(const std::string &name1, const std::string &name2);
    void reset(const std::string &name);
    double acc(const std::string &name);

    void report() const;

    const maptype &table() const { return lut; }
  };

extern wallTimerSet wallTimers;

void tstack_push(const std::string &name);
void tstack_pop(const std::string &name);
void tstack_pop();
void tstack_replace(const std::string &name1, const std::string &name2);
void tstack_replace(const std::string &name);
void tstack_report(const std::string &stem);

#endif
