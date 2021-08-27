/*
 *  This file is part of the MR utility library.
 *
 *  This code is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This code is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this code; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/** \file ducc0/math/gl_integrator.h
 *  Functionality for Gauss-Legendre quadrature
 *
 *  \copyright Copyright (C) 2019-2021 Max-Planck-Society, Ignace Bogaert
 *  \author Martin Reinecke
 */

#ifndef DUCC0_GL_INTEGRATOR_H
#define DUCC0_GL_INTEGRATOR_H

#include <cmath>
#include <array>
#include <cstddef>
#include <utility>
#include <vector>
#include "ducc0/math/constants.h"
#include "ducc0/infra/error_handling.h"

namespace ducc0 {

namespace detail_gl_integrator {

using namespace std;

template<typename T> inline T one_minus_x2 (T x)
  { if (x<0) x=-x; return (x>T(0.1)) ? (T(1)+x)*(T(1)-x) : T(1)-x*x; }

pair<double, double> calc_gl_iterative(size_t n, size_t i)
  {
  using Tfloat = long double;
  constexpr Tfloat eps=Tfloat(3e-14L);
  const Tfloat dn = Tfloat(n);
  const Tfloat t0 = Tfloat(1) - (1-Tfloat(1)/dn) / (Tfloat(8)*dn*dn);
  const Tfloat t1 = Tfloat(1)/(Tfloat(4)*dn+Tfloat(2));
  Tfloat x0 = cos(double(pi * ((i<<2)-1) * t1)) * t0;

  bool dobreak=false;
  size_t j=0;
  Tfloat dpdx;
  while(1)
    {
    Tfloat P_1 = 1;
    Tfloat P0 = x0;
    Tfloat dx, x1;

    for (size_t k=2; k<=n; k++)
      {
      Tfloat P_2 = P_1;
      P_1 = P0;
      P0 = x0*P_1 + (Tfloat(k)-Tfloat(1))/Tfloat(k) * (x0*P_1-P_2);
      }

    dpdx = (P_1 - x0*P0) * n / one_minus_x2(x0);

    /* Newton step */
    x1 = x0 - P0/dpdx;
    dx = x0-x1;
    x0 = x1;
    if (dobreak) break;

    if ((dx<0?-dx:dx)<=eps) dobreak=1;
    MR_assert(++j<100, "convergence problem");
    }

  return make_pair(double(x0), double(2./(one_minus_x2(x0)*dpdx*dpdx)));
  }


// The next three functions are modified versions of the FastGL code
// by Ignace Bogaert. The code is available at
// https://sourceforge.net/projects/fastgausslegendrequadrature/ 
// A paper describing the algorithms is available at
// https://epubs.siam.org/doi/pdf/10.1137/140954969

// This function computes the kth zero of the BesselJ(0,x)
double besseljzero(int k)
  {
  constexpr static array<double,20> JZ
    {2.40482555769577276862163187933,  5.52007811028631064959660411281,
     8.65372791291101221695419871266,  11.7915344390142816137430449119,
     14.9309177084877859477625939974, 18.0710639679109225431478829756,
     21.2116366298792589590783933505, 24.3524715307493027370579447632,
     27.4934791320402547958772882346, 30.6346064684319751175495789269,
     33.7758202135735686842385463467, 36.9170983536640439797694930633,
     40.0584257646282392947993073740, 43.1997917131767303575240727287,
     46.3411883716618140186857888791, 49.4826098973978171736027615332,
     52.6240518411149960292512853804, 55.7655107550199793116834927735,
     58.9069839260809421328344066346, 62.0484691902271698828525002646};

  if (k<=20) return JZ[k-1];

  double z = pi*(k-0.25);
  double r = 1.0/z;
  double r2 = r*r;
  return z
       + r*(0.125
       +r2*(-0.807291666666666666666666666667e-1
       +r2*(0.246028645833333333333333333333
       +r2*(-1.82443876720610119047619047619
       +r2*(25.3364147973439050099206349206
       +r2*(-567.644412135183381139802038240
       +r2*(18690.4765282320653831636345064
       +r2*(-8.49353580299148769921876983660e5
       +r2*5.09225462402226769498681286758e7))))))));
  }

// This function computes the square of BesselJ(1, BesselZero(0,k))
double besselj1squared(int k)
  {
  constexpr static array<double,21> J1
    {0.269514123941916926139021992911 , 0.115780138582203695807812836182,
     0.0736863511364082151406476811985, 0.0540375731981162820417749182758,
     0.0426614290172430912655106063495, 0.0352421034909961013587473033648,
     0.0300210701030546726750888157688, 0.0261473914953080885904584675399,
     0.0231591218246913922652676382178, 0.0207838291222678576039808057297,
     0.0188504506693176678161056800214, 0.0172461575696650082995240053542,
     0.0158935181059235978027065594287, 0.0147376260964721895895742982592,
     0.0137384651453871179182880484134, 0.0128661817376151328791406637228,
     0.0120980515486267975471075438497, 0.0114164712244916085168627222986,
     0.0108075927911802040115547286830, 0.0102603729262807628110423992790,
     0.00976589713979105054059846736696};

  if (k<=21) return J1[k-1];

  double x = 1.0/(k-0.25);
  double x2 = x*x;
  return x * (0.202642367284675542887758926420
       + x2*x2*(-0.303380429711290253026202643516e-3
       + x2*(0.198924364245969295201137972743e-3
       + x2*(-0.228969902772111653038747229723e-3
       + x2*(0.433710719130746277915572905025e-3
       + x2*(-0.123632349727175414724737657367e-2
       + x2*(0.496101423268883102872271417616e-2
       + x2*(-0.266837393702323757700998557826e-1
       + x2*.185395398206345628711318848386))))))));
  }

// Compute a node-weight pair, with k limited to half the range
pair<double,double> calc_gl_bogaert(size_t n, size_t k0)
  {
  size_t k = ((2*k0-1)<=n) ? k0 : n-k0+1;
  // First get the Bessel zero
  double w = 1.0/(n+0.5);
  double nu = besseljzero(k);
  double theta = w*nu;
  double x = theta*theta;

  // Get the asymptotic BesselJ(1,nu) squared
  double B = besselj1squared(k);

  // Get the Chebyshev interpolants for the nodes...
  double SF1T = (((((-1.29052996274280508473467968379e-12*x
                     +2.40724685864330121825976175184e-10)*x
                     -3.13148654635992041468855740012e-8)*x
                     +0.275573168962061235623801563453e-5)*x
                     -0.148809523713909147898955880165e-3)*x
                     +0.416666666665193394525296923981e-2)*x
                     -0.416666666666662959639712457549e-1;
  double SF2T = (((((+2.20639421781871003734786884322e-9*x
                     -7.53036771373769326811030753538e-8)*x
                     +0.161969259453836261731700382098e-5)*x
                     -0.253300326008232025914059965302e-4)*x
                     +0.282116886057560434805998583817e-3)*x
                     -0.209022248387852902722635654229e-2)*x
                     +0.815972221772932265640401128517e-2;
  double SF3T = (((((-2.97058225375526229899781956673e-8*x
                     +5.55845330223796209655886325712e-7)*x
                     -0.567797841356833081642185432056e-5)*x
                     +0.418498100329504574443885193835e-4)*x
                     -0.251395293283965914823026348764e-3)*x
                     +0.128654198542845137196151147483e-2)*x
                     -0.416012165620204364833694266818e-2;

  // ...and for the weights
  double WSF1T = ((((((((-2.20902861044616638398573427475e-14*x
                         +2.30365726860377376873232578871e-12)*x
                         -1.75257700735423807659851042318e-10)*x
                         +1.03756066927916795821098009353e-8)*x
                         -4.63968647553221331251529631098e-7)*x
                         +0.149644593625028648361395938176e-4)*x
                         -0.326278659594412170300449074873e-3)*x
                         +0.436507936507598105249726413120e-2)*x
                         -0.305555555555553028279487898503e-1)*x
                         +0.833333333333333302184063103900e-1;
  double WSF2T = (((((((+3.63117412152654783455929483029e-12*x
                        +7.67643545069893130779501844323e-11)*x
                        -7.12912857233642220650643150625e-9)*x
                        +2.11483880685947151466370130277e-7)*x
                        -0.381817918680045468483009307090e-5)*x
                        +0.465969530694968391417927388162e-4)*x
                        -0.407297185611335764191683161117e-3)*x
                        +0.268959435694729660779984493795e-2)*x
                        -0.111111111111214923138249347172e-1;
  double WSF3T = (((((((+2.01826791256703301806643264922e-9*x
                        -4.38647122520206649251063212545e-8)*x
                        +5.08898347288671653137451093208e-7)*x
                        -0.397933316519135275712977531366e-5)*x
                        +0.200559326396458326778521795392e-4)*x
                        -0.422888059282921161626339411388e-4)*x
                        -0.105646050254076140548678457002e-3)*x
                        -0.947969308958577323145923317955e-4)*x
                        +0.656966489926484797412985260842e-2;

  // Then refine with the paper expansions
  double NuoSin = nu/sin(theta);
  double BNuoSin = B*NuoSin;
  double WInvSinc = w*w*NuoSin;
  double WIS2 = WInvSinc*WInvSinc;

  // Finally compute the node and the weight
  theta = w*(nu + theta * WInvSinc * (SF1T + WIS2*(SF2T + WIS2*SF3T)));
  double Deno = BNuoSin + BNuoSin * WIS2*(WSF1T + WIS2*(WSF2T + WIS2*WSF3T));
  double weight = (2.0*w)/Deno;
  return make_pair((k==k0) ? cos(theta) : -cos(theta), weight);
  }

pair<double, double> calc_gl(size_t n, size_t k)
  {
  MR_assert(n>=k, "k must not be greater than n");
  MR_assert(k>0, "k must be positive");
  return (n<=100) ? calc_gl_iterative(n,k) : calc_gl_bogaert(n,k);
  }

/// Class for computing Gauss-Lgendre abscissas, weights and intgrals
class GL_Integrator
  {
  private:
    size_t n_;
    vector<double> x, w;

  public:
    /// Creates an integrator for \a n abscissas
    /** \note The \a nthreads parameter is obsolescent and ignored. */
    GL_Integrator(size_t n, size_t /*nthreads*/=1)
      : n_(n)
      {
      MR_assert(n>=1, "number of points must be at least 1");
      size_t m = (n+1)>>1;
      x.resize(m);
      w.resize(m);
      for (size_t i=0; i<m; ++i)
        {
        auto tmp = calc_gl(n, m-i);
        x[i] = tmp.first;
        w[i] = tmp.second;
        }
      }

    /// Returns the approximated integral of \a func in [-1;1]
    template<typename Func> auto integrate(Func f) -> decltype(f(0.))
      {
      using T = decltype(f(0.));
      T res=0;
      size_t istart=0;
      if (n_&1)
        {
        res = f(x[0])*w[0];
        istart=1;
        }
      for (size_t i=istart; i<x.size(); ++i)
        res += (f(x[i])+f(-x[i]))*w[i];
      return res;
      }

    /// Returns the approximated integral of \a func in [-1;1], where \a func
    /// is symmetric with respect to x=0.
    template<typename Func> auto integrateSymmetric(Func f) -> decltype(f(0.))
      {
      using T = decltype(f(0.));
      T res=f(x[0])*w[0];
      if (n_&1) res *= 0.5;
      for (size_t i=1; i<x.size(); ++i)
        res += f(x[i])*w[i];
      return res*2;
      }

    /// Returns the Gauss-Legendre abscissas.
    vector<double> coords() const
      {
      vector<double> res(n_);
      for (size_t i=0; i<x.size(); ++i)
        {
        res[i]=-x[x.size()-1-i];
        res[n_-1-i] = x[x.size()-1-i];
        }
      return res;
      }
    /// Returns the non-negative Gauss-Legendre abscissas.
    const vector<double> &coordsSymmetric() const
      { return x; }

    /// Returns the Gauss-Legendre weights.
    vector<double> weights() const
      {
      vector<double> res(n_);
      for (size_t i=0; i<w.size(); ++i)
        res[i]=res[n_-1-i]=w[w.size()-1-i];
      return res;
      }
    /// Returns the Gauss-Legendre weights for the non-negative abscissas,
    /// with an additional factor of 2 for positive abscissas.
    vector<double> weightsSymmetric() const
      {
      auto res = w;
      if (n_&1) res[0]*=0.5;
      for (auto &v:res) v*=2;
      return res;
      }
  };

}

using detail_gl_integrator::GL_Integrator;

}

#endif
