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

#ifndef BINOMIAL_H
#define BINOMIAL_H

#include <assert.h>
#include <vector>

using std::vector;

/**
 * A template class facilitating fast computation and retrieval of binomials.
 * The binomials are computed only once at first call of Get function according
 * to Pascal's Triangle.
 */
template <class T> class Binomial {
public:
  typedef vector<T> VectorT;
  typedef vector<VectorT> VVectorT;

  /** Retrieves the binomial _i "over" _j */
  static T Get(int _i, int _j);
  /** Sets the maximal value of upper binomial param to _max */
  static void SetMax(int _max);
  /** Gets the maximal value of upper binomial param */
  static int GetMax();

private:
  /** Computed Pascal's Triangle */
  static void ComputePascalsTriangle();

  static VVectorT pascalsTriangle_;
  static int max_;
};

template <class T> inline void Binomial<T>::SetMax(int _max) {
  max_ = _max;
  ComputePascalsTriangle();
}

template <class T> inline int Binomial<T>::GetMax() { return max_; }

template <class T> typename Binomial<T>::VVectorT Binomial<T>::pascalsTriangle_;

template <class T> int Binomial<T>::max_ = 60;

// Moved from Binomial.cpp to avoid templated instantation
template <class T> void Binomial<T>::ComputePascalsTriangle() {
  // allocate storage for the pascal triangle of size determined by max_
  pascalsTriangle_.resize(max_ + 1);

  for (int i = 0; i < max_ + 1; ++i) {
    pascalsTriangle_[i].resize(max_ + 1 - i);
    for (int j = 0; j < max_ + 1 - i; ++j) {
      // the values are ones on the edges of the triangle
      if (!i || !j) {
        pascalsTriangle_[i][j] = (T)1;
      }
      // use the familiar addition to generate values on lower levels
      else {
        pascalsTriangle_[i][j] =
            pascalsTriangle_[i][j - 1] + pascalsTriangle_[i - 1][j];
      }
    }
  }
}

template <class T> T Binomial<T>::Get(int _i, int _j) {
  // the values are computed only first time this function is called
  if (!pascalsTriangle_.size()) {
    ComputePascalsTriangle();
  }

  assert(_i >= 0 && _j >= 0 && _i >= _j);

  return pascalsTriangle_[_j][_i - _j];
}

#endif
