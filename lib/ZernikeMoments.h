/*

                          3D Zernike Moments
    Copyright (C) 2003 by Computer Graphics Group, University of Bonn
           http://www.cg.cs.uni-bonn.de/project-pages/3dsearch/

Code by Marcin Novotni:     marcin@cs.uni-bonn.de

for more information, see the paper:

@inproceedings{novotni-2003-3d,
    author = {M. Novotni and R. Klein},
    title = {3{D} {Z}ernike Descriptors for Content Based Shape Retrieval},
    booktitle = {The 8th ACM Symposium on Solid Modeling and Applications},
    pages = {216--225},
    year = {2003},
    month = {June},
    institution = {Universit\"{a}t Bonn},
    conference = {The 8th ACM Symposium on Solid Modeling and Applications, June
16-20, Seattle, WA}
}
 *---------------------------------------------------------------------------*
 *                                                                           *
 *                                License                                    *
 *                                                                           *
 *  This library is free software; you can redistribute it and/or modify it  *
 *  under the terms of the GNU Library General Public License as published   *
 *  by the Free Software Foundation, version 2.                              *
 *                                                                           *
 *  This library is distributed in the hope that it will be useful, but      *
 *  WITHOUT ANY WARRANTY; without even the implied warranty of               *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU        *
 *  Library General Public License for more details.                         *
 *                                                                           *
 *  You should have received a copy of the GNU Library General Public        *
 *  License along with this library; if not, write to the Free Software      *
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                *
 *                                                                           *
\*===========================================================================*/

#pragma once

#define PI 3.141592653589793

#include <complex>
#include <cstdlib>
#include <ios>
#include <omp.h>
#include <set>

#include "Binomial.h"
#include "Factorial.h"
#include "ScaledGeometricMoments.h"

void printProgressBar(int progress, int total, int barWidth = 60) {
  float ratio = static_cast<float>(progress) / total;
  int filledWidth = static_cast<int>(barWidth * ratio);

  std::cout << "\r[";
  for (int i = 0; i < barWidth; ++i) {
    if (i < filledWidth)
      std::cout << "#";
    else
      std::cout << "-";
  }
  std::cout << "] " << static_cast<int>(ratio * 100.0) << "%";
  std::cout.flush();
}

/**
 * Struct representing a complex coefficient of a moment of order (p_,q_,r_).
 */
template <class T> struct ComplexCoeff {
  typedef std::complex<T> ValueT;

  // Constructor with scalar args
  ComplexCoeff(int p, int q, int r, const ValueT &value)
      : p_(p), q_(q), r_(r), value_(value) {}

  // Copy constructor
  ComplexCoeff(const ComplexCoeff<T> &cc)
      : p_(cc.p_), q_(cc.q_), r_(cc.r_), value_(cc.value_) {}

  // Default constructor
  ComplexCoeff() : p_(0), q_(0), r_(0) {}

  int p_, q_, r_;
  ValueT value_;
};

/**
 * Class representing the Zernike moments
 */
template <class VoxelT, class MomentT> class ZernikeMoments {
public:
  typedef MomentT T;
  typedef vector<T> T1D;   // vector of scalar type
  typedef vector<T1D> T2D; // 2D array of scalar type
  typedef vector<T2D> T3D; // 3D array of scalar type
  typedef vector<T3D> T4D; // 3D array of scalar type

  typedef std::complex<T> ComplexT; // complex type
  typedef vector<vector<vector<ComplexT>>>
      ComplexT3D; // 3D array of complex type

  typedef ComplexCoeff<T> ComplexCoeffT;
  typedef vector<vector<vector<vector<ComplexCoeffT>>>> ComplexCoeffT4D;

public:
  ZernikeMoments() : order_(0) {}
  ZernikeMoments(int _order, ScaledGeometricalMoments<VoxelT, MomentT> &_gm) {
    Init(_order, _gm);
  }

  /**
   * Computes all coefficients that are input data independent
   */
  void Init(int _order, ScaledGeometricalMoments<VoxelT, MomentT> &_gm) {
    gm_ = _gm;
    order_ = _order;

    ComputeCs();
    ComputeQs();
    ComputeGCoefficients();
  }

  /**
   * Computes the Zernike moments for current object.
   */
  void Compute() {
    // geometrical moments have to be computed first
    if (!order_) {
      std::cerr
          << "ZernikeMoments<VoxelT,MomentT>::ComputeZernikeMoments (): attempting to \
                       compute Zernike moments without setting valid geometrical \
                       moments first. \nExiting...\n";
      exit(-1);
    }

    /*
     indexing:
       n goes 0..order_
       l goes 0..n so that n-l is even
       m goes -l..l
    */

    T nullMoment;
    zernikeMoments_.resize(order_ + 1);
    for (int n = 0; n <= order_; ++n) {
      zernikeMoments_[n].resize(n / 2 + 1);

      int l0 = n % 2, li = 0;
      for (int l = l0; l <= n; ++li, l += 2) {
        zernikeMoments_[n][li].resize(l + 1);
        for (int m = 0; m <= l; ++m) {
          // Zernike moment of according indices [nlm]
          ComplexT zm((T)0, (T)0);

          int nCoeffs = (int)gCoeffs_[n][li][m].size();
          for (int i = 0; i < nCoeffs; ++i) {
            ComplexCoeffT cc = gCoeffs_[n][li][m][i];
            zm += std::conj(cc.value_) * gm_.GetMoment(cc.p_, cc.q_, cc.r_);
          }

          zm *= (T)(3.0 / (4.0 * PI));
          if (n == 0 && l == 0 && m == 0) {
            // FixMe Unused variable! What is it for?
            nullMoment = zm.real();
          }
          zernikeMoments_[n][li][m] = zm;
        }
      }
    }
  }

  /**
   *  Access Zernike Moments.
   *
   * @param _n Order along x.
   * @param _l Order along y.
   * @param _m Order along z.
   * @return Complex Moment.
   */
  inline ComplexT GetMoment(int _n, int _l, int _m) {
    if (_m >= 0) {
      return zernikeMoments_[_n][_l / 2][_m];
    } else {
      T sign;
      if (_m % 2) {
        sign = (T)(-1);
      } else {
        sign = (T)1;
      }
      return sign * std::conj(zernikeMoments_[_n][_l / 2][abs(_m)]);
    }
  }
  /**
   *  Reconstructs the original object from the computed moments.
   *
   * @param _grid Complex output grid.
   * @param _xCOG Scaled center of gravity.
   * @param _yCOG Scaled center of gravity.
   * @param _zCOG Scaled center of gravity.
   * @param _minN Min value for n freq index.
   * @param _maxN Max value for n freq index.
   * @param _minL Max value for l freq index.
   * @param _maxL Max value for l freq index.
   */
  void Reconstruct(ComplexT3D &_grid, T _xCOG, T _yCOG, T _zCOG, T _scale,
                   int _minN, int _maxN, int _minL, int _maxL) {
    if (_maxN == 100)
      _maxN = order_;

    int dimX = _grid.size();
    int dimY = _grid[0].size();
    int dimZ = _grid[0][0].size();

    // translation
    T vx = _xCOG;
    T vy = _yCOG;
    T vz = _zCOG;

    // Unit sphere constraints
    T dxy2, dxyz2;
    std::vector<T> dx(dimX), dx2(dimX);
    std::vector<T> dy(dimY), dy2(dimY);
    std::vector<T> dz(dimZ), dz2(dimZ);

    for (int i = 0; i < dimX; ++i) {
      dx[i] = ((T)i - vx) * _scale;
      dx2[i] = dx[i] * dx[i];
    }
    for (int i = 0; i < dimY; ++i) {
      dy[i] = ((T)i - vy) * _scale;
      dy2[i] = dy[i] * dy[i];
    }
    for (int i = 0; i < dimZ; ++i) {
      dz[i] = ((T)i - vz) * _scale;
      dz2[i] = dz[i] * dz[i];
    }

    long total_iterations = dimX;
    long progress = 0;
// origin is at the grid center, all voxels are projected onto the unit sphere
#pragma omp parallel for schedule(dynamic)
    for (int x = 0; x < dimX; ++x) {

#pragma omp critical
      { printProgressBar(progress, total_iterations); }

      for (int y = 0; y < dimY; ++y) {
        dxy2 = dx2[y] + dy2[x];
        for (int z = 0; z < dimZ; ++z) {
          dxyz2 = dxy2 + dz2[z];

          if (dxyz2 > 1.0)
            continue;

          // function value at point
          ComplexT fVal = (0, 0);
          for (int n = _minN; n <= _maxN; ++n) {
            int maxK = n / 2;
            for (int k = 0; k <= maxK; ++k) {
              for (int nu = 0; nu <= k; ++nu) {
                // check whether l is within bounds
                int l = n - 2 * k;
                if (l < _minL || l > _maxL)
                  continue;

                for (int m = -l; m <= l; ++m) {
                  // zernike polynomial evaluated at point
                  ComplexT zp(0, 0);

                  int absM = std::abs(m);

                  int nCoeffs = gCoeffs_[n][l / 2][absM].size();
                  for (int i = 0; i < nCoeffs; ++i) {
                    ComplexCoeffT cc = gCoeffs_[n][l / 2][absM][i];
                    ComplexT cvalue = cc.value_;

                    // conjugate if m negative and take care of sign
                    if (m < 0) {
                      cvalue = std::conj(cvalue);
                      if (m % 2)
                        cvalue *= (T)(-1);
                    }

                    zp += cvalue * std::pow(dx[x], (T)cc.p_) *
                          std::pow(dy[y], (T)cc.q_) * std::pow(dz[z], (T)cc.r_);
                  }

                  fVal += zp * GetMoment(n, l, m);
                }
              }
            }
          }
          _grid[x][y][z] = fVal;
        }
      }
      progress++;
    }

    // NormalizeGridValues (_grid);
    printProgressBar(total_iterations, total_iterations);
    std::cout << std::endl;
  }

  void NormalizeGridValues(ComplexT3D &_grid);
  void CheckOrthonormality(int _n1, int _l1, int _m1, int _n2, int _l2,
                           int _m2);

private:
  void ComputeCs();
  void ComputeQs();
  void ComputeGCoefficients();
  T EvalMonomialIntegral(int _p, int _q, int _r, int _dim);

private:
  ComplexCoeffT4D gCoeffs_;   // coefficients of the geometric moments
  ComplexT3D zernikeMoments_; // nomen est omen
  T3D qs_; // q coefficients (radial polynomial normalization)
  T2D cs_; // c coefficients (harmonic polynomial normalization)

  ScaledGeometricalMoments<VoxelT, MomentT> gm_;
  int order_; // := max{n} according to indexing of Zernike polynomials
};

/**
 * Computes all the normalizing factors $c_l^m$ for harmonic polynomials
 * Indexing
 *  l goes from 0 to n
 *  m goes from -l to l, in fact from 0 to l, since c(l,-m) = c (l,m)
 */
template <class VoxelT, class MomentT>
void ZernikeMoments<VoxelT, MomentT>::ComputeCs() {
  Factorial<T> factorial;
  cs_.resize(order_ + 1);
  for (int l = 0; l <= order_; ++l) {
    cs_[l].resize(l + 1);
    for (int m = 0; m <= l; ++m) {
      T n_sqrt = ((T)2 * l + (T)1) * factorial.Get(l + 1, l + m);
      T d_sqrt = factorial.Get(l - m + 1, l);
      cs_[l][m] = sqrt(n_sqrt / d_sqrt);
    }
  }
}

/**
 * Computes all coefficients q for orthonormalization of radial polynomials
 * in Zernike polynomials.
 * Indexing
 *  n goes 0..order_
 *  l goes 0..n, so that n-l is even
 *  mu goes 0..(n-l)/2
 */
template <class VoxelT, class MomentT>
void ZernikeMoments<VoxelT, MomentT>::ComputeQs() {
  Binomial<T> binomial;
  qs_.resize(order_ + 1); // there is order_ + 1 n's
  for (int n = 0; n <= order_; ++n) {
    qs_[n].resize(n / 2 + 1); // there is floor(n/2) + 1 l's

    int l0 = n % 2;
    for (int l = l0; l <= n; l += 2) {
      int k = (n - l) / 2;

      qs_[n][l / 2].resize(k + 1); // there is k+1 mu's
      for (int mu = 0; mu <= k; ++mu) {
        T nom = binomial.Get(2 * k, k) * // nominator of straight part
                binomial.Get(k, mu) * binomial.Get(2 * (k + l + mu) + 1, 2 * k);

        if ((k + mu) % 2) {
          nom *= (T)(-1);
        }

        // denominator of straight part
        T den = std::pow((T)2, (T)(2 * k)) * binomial.Get(k + l + mu, k);
        T n_sqrt = (T)(2 * l + 4 * k + 3); // nominator of sqrt part
        T d_sqrt = (T)3;                   // denominator of sqrt part

        qs_[n][l / 2][mu] = nom / den * sqrt(n_sqrt / d_sqrt);
      }
    }
  }
}

/**
 * Computes the coefficients of geometrical moments in linear combinations
 * yielding the Zernike moments for each applicable [n,l,m] for n<=order_.
 * For each such combination the coefficients are stored with according
 * geom. moment order (see ComplexCoeff).
 */
template <class VoxelT, class MomentT>
void ZernikeMoments<VoxelT, MomentT>::ComputeGCoefficients() {
  int countCoeffs = 0;
  gCoeffs_.resize(order_ + 1);

  Binomial<T> binomial;
  for (int n = 0; n <= order_; ++n) {
    gCoeffs_[n].resize(n / 2 + 1);
    int li = 0, l0 = n % 2;
    for (int l = l0; l <= n; ++li, l += 2) {
      gCoeffs_[n][li].resize(l + 1);
      for (int m = 0; m <= l; ++m) {
        T w = cs_[l][m] / std::pow((T)2, (T)m);

        int k = (n - l) / 2;
        for (int nu = 0; nu <= k; ++nu) {
          T w_Nu = w * qs_[n][li][nu];
          for (int alpha = 0; alpha <= nu; ++alpha) {
            T w_NuA = w_Nu * binomial.Get(nu, alpha);
            for (int beta = 0; beta <= nu - alpha; ++beta) {
              T w_NuAB = w_NuA * binomial.Get(nu - alpha, beta);
              for (int p = 0; p <= m; ++p) {
                T w_NuABP = w_NuAB * binomial.Get(m, p);
                for (int mu = 0; mu <= (l - m) / 2; ++mu) {
                  T w_NuABPMu = w_NuABP * binomial.Get(l, mu) *
                                binomial.Get(l - mu, m + mu) /
                                (T)pow(2.0, (double)(2 * mu));
                  for (int q = 0; q <= mu; ++q) {
                    // the absolute value of the coefficient
                    T w_NuABPMuQ = w_NuABPMu * binomial.Get(mu, q);

                    // the sign
                    if ((m - p + mu) % 2) {
                      w_NuABPMuQ *= T(-1);
                    }

                    // * i^p
                    int rest = p % 4;
                    ComplexT c;
                    switch (rest) {
                    case 0:
                      c = ComplexT(w_NuABPMuQ, (T)0);
                      break;
                    case 1:
                      c = ComplexT((T)0, w_NuABPMuQ);
                      break;
                    case 2:
                      c = ComplexT((T)(-1) * w_NuABPMuQ, (T)0);
                      break;
                    case 3:
                      c = ComplexT((T)0, (T)(-1) * w_NuABPMuQ);
                      break;
                    }

                    // determination of the order of according moment
                    int z_i = l - m + 2 * (nu - alpha - beta - mu);
                    int y_i = 2 * (mu - q + beta) + m - p;
                    int x_i = 2 * q + p + 2 * alpha;

                    ComplexCoeffT cc(x_i, y_i, z_i, c);
                    gCoeffs_[n][li][m].push_back(cc);
                    countCoeffs++;

                  } // q
                } // mu
              } // p
            } // beta
          } // alpha
        } // nu
      } // m
    } // l
  } // n
}

template <class VoxelT, class MomentT>
void ZernikeMoments<VoxelT, MomentT>::NormalizeGridValues(ComplexT3D &_grid) {
  int xD = _grid.size();
  int yD = _grid[0].size();
  int zD = _grid[0][0].size();

  T max = (T)0;
  for (int k = 0; k < zD; ++k) {
    for (int j = 0; j < yD; ++j) {
      for (int i = 0; i < xD; ++i) {
        if (_grid[i][j][k].real() > max) {
          max = _grid[i][j][k].real();
        }
      }
    }
  }

  std::cout << "\nMaximal value in grid: " << max << "\n";
  T invMax = (T)1 / max;
  for (int k = 0; k < zD; ++k) {
    for (int j = 0; j < yD; ++j) {
      for (int i = 0; i < xD; ++i) {
        _grid[i][j][k] *= invMax;
      }
    }
  }
}

template <class VoxelT, class MomentT>
void ZernikeMoments<VoxelT, MomentT>::CheckOrthonormality(int _n1, int _l1,
                                                          int _m1, int _n2,
                                                          int _l2, int _m2) {
  int li1 = _l1 / 2;
  int li2 = _l2 / 2;
  int dim = 64;

  // the total sum of the scalar product
  ComplexT sum((T)0, (T)0);

  int nCoeffs1 = (int)gCoeffs_[_n1][li1][_m1].size();
  int nCoeffs2 = (int)gCoeffs_[_n2][li2][_m2].size();

  for (int i = 0; i < nCoeffs1; ++i) {
    ComplexCoeffT cc1 = gCoeffs_[_n1][li1][_m1][i];
    for (int j = 0; j < nCoeffs2; ++j) {
      ComplexCoeffT cc2 = gCoeffs_[_n2][li2][_m2][j];

      T temp = (T)0;

      int p = cc1.p_ + cc2.p_;
      int q = cc1.q_ + cc2.q_;
      int r = cc1.r_ + cc2.r_;

      sum += cc1.value_ * std::conj(cc2.value_) *
             EvalMonomialIntegral(p, q, r, dim);
    }
  }

  std::cout << "\nInner product of [" << _n1 << "," << _l1 << "," << _m1 << "]";
  std::cout << " and [" << _n2 << "," << _l2 << "," << _m2 << "]: ";
  std::cout << sum << "\n\n";
}

/**
 * Evaluates the integral of a monomial x^p*y^q*z^r within the unit sphere
 * Attention : a very stupid implementation, thus it's accordingly very slow
 */
template <class VoxelT, class MomentT>
MomentT ZernikeMoments<VoxelT, MomentT>::EvalMonomialIntegral(int _p, int _q,
                                                              int _r,
                                                              int _dim) {
  T radius = (T)(_dim - 1) / (T)2;
  T scale = std::pow((T)1 / radius, 3);
  T center = (T)(_dim - 1) / (T)2;

  T result = (T)0;
  T point[3];

  for (int x = 0; x < _dim; ++x) {
    point[0] = ((T)x - center) / radius;
    for (int y = 0; y < _dim; ++y) {
      point[1] = ((T)y - center) / radius;
      for (int z = 0; z < _dim; ++z) {
        point[2] = ((T)z - center) / radius;

        if (point[0] * point[0] + point[1] * point[1] + point[2] * point[2] >
            (T)1) {
          continue;
        }

        result += std::pow(point[0], (T)_p) * std::pow(point[1], (T)_q) *
                  std::pow(point[2], (T)_r);
      }
    }
  }

  result *= (T)(3.0 / (4.0 * PI)) * scale;
  return result;
}
