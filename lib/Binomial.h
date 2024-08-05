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

#include <assert.h>
#include <vector>

/**
 * Class facilitating fast computation and retrieval of binomials using Pascal's
 * triangle.
 */
template <class T> class Binomial {
public:
  typedef std::vector<T> VectorT;
  typedef std::vector<VectorT> VVectorT;

public:
  /** Retrieves the binomial _i "over" _j */
  static T Get(int _i, int _j) {
    if (!pascalsTriangle_.size()) {
      ComputePascalsTriangle();
    }
    assert(_i >= 0 && _j >= 0 && _i >= _j);
    return pascalsTriangle_[_j][_i - _j];
  }

  /** Sets the maximal value of upper binomial param to _max */
  static inline void SetMax(int _max) {
    max_ = _max;
    ComputePascalsTriangle();
  };
  /** Gets the maximal value of upper binomial param */
  static inline int GetMax() { return max_; };

private:
  /** Computed Pascal's Triangle */
  static void ComputePascalsTriangle() {
    pascalsTriangle_.resize(max_ + 1);

    for (int i = 0; i < max_ + 1; ++i) {
      pascalsTriangle_[i].resize(max_ + 1 - i);
      for (int j = 0; j < max_ + 1 - i; ++j) {
        // the values are ones on the edges of the triangle
        if (!i || !j) {
          pascalsTriangle_[i][j] = (T)1;
        } else {
          pascalsTriangle_[i][j] =
              pascalsTriangle_[i][j - 1] + pascalsTriangle_[i - 1][j];
        }
      }
    }
  }

private:
  static VVectorT pascalsTriangle_;
  static int max_;
};

template <class T> typename Binomial<T>::VVectorT Binomial<T>::pascalsTriangle_;

template <class T> int Binomial<T>::max_ = 60;
