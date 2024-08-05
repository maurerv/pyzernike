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
#include <iostream>
#include <vector>

using std::vector;

/**
 * Class for precomputation and subsequent retrieval of factorials of an integer
 * number.
 */
template <class T> class Factorial {
public:
  /** Gets the factorial of _i */
  static inline T Get(int _i) {
    assert(_i >= 0 && _i <= max_);

    if (!factorials_.size()) {
      ComputeFactorials();
    }

    // 0! = 1
    if (!_i) {
      return 1;
    }
    return factorials_[_i - 1];
  }

  /** Gets _i*(_i+1)*...*(_j-1)*_j */
  static inline T Get(int _i, int _j) {
    T result = (T)1;

    for (int i = _j; i >= _i; --i) {
      result *= i;
    }
    return result;
  }

  /** Sets the maximal stored factorial value to _max */
  static inline void SetMax(int _max) {
    assert(_max >= (T)0);
    max_ = _max;
    if (max_ <= _max) {
      ComputeFactorials();
    }
  }
  /** Gets the maximal stored factorial value */
  static inline int GetMax() { return max_; };

private:
  /** Computes factorials of numbers [1..max] */
  static inline void ComputeFactorials() {
    factorials_.resize(max_);
    factorials_[0] = (T)1;
    for (int i = 1; i < max_; ++i) {
      factorials_[i] = factorials_[i - 1] * (T)(i + 1);
    }
  }

  static int max_;
  static vector<T> factorials_;
  static vector<vector<T>> subFactorials_;
};

// The obligatory initialization of static attributes
template <class T> int Factorial<T>::max_ = 19;

template <class T> vector<T> Factorial<T>::factorials_;
