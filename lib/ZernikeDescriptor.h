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

#include <fstream>
#include <iostream>
#include <vector>

#include "ScaledGeometricMoments.h"
#include "ZernikeMoments.h"

/**
 * Interface for dealing with 3D Zernike Descriptors.
 */
template <class T, class TIn> class ZernikeDescriptor {
public:
  typedef std::complex<T> ComplexT;
  typedef vector<vector<vector<ComplexT>>> ComplexT3D;

  typedef vector<T> T1D;
  typedef vector<T1D> T2D;

  typedef ZernikeMoments<T, T> ZernikeMomentsT;
  typedef ScaledGeometricalMoments<T, T> ScaledGeometricalMomentsT;

public:
  /**
   * Initialize ZernikeDescriptor from file.
   *
   * @param rawName Path to binary file containing cubic grid.
   * @param order Maximal order of the moments.
   * @return ZernikeDescriptor instance.
   */
  ZernikeDescriptor(const char *rawName, int order)
      : ZernikeDescriptor(ReadGrid(rawName, dim_), dim_, order) {}

  /**
   * Initialize ZernikeDescriptor from grid.
   *
   * @param voxels Cubic voxel grid.
   * @order dim Dimension of the grid.
   * @order order Maximal order of the moments.
   * @return ZernikeDescriptor instance.
   */
  ZernikeDescriptor(T *voxels, int dim, int order)
      : voxels_(voxels), dim_(dim), order_(order) {
    ComputeNormalization();
    NormalizeGrid();
    ComputeMoments();
    ComputeInvariants();
  }

  /**
   *  Reconstructs the original object from the computed moments.
   *
   * @param _grid Complex output grid.
   * @param _minN Min value for n freq index.
   * @param _maxN Max value for n freq index.
   * @param _minL Max value for l freq index.
   * @param _maxL Max value for l freq index.
   *
   */
  void Reconstruct(ComplexT3D &_grid, int _minN = 0, int _maxN = 100,
                   int _minL = 0, int _maxL = 100) {
    // the scaling between the reconstruction and original grid
    T fac = (T)(_grid.size()) / (T)dim_;
    zm_.Reconstruct(_grid, center[0] * fac, center[1] * fac, center[2] * fac,
                    scale_ / fac, _minN, _maxN, _minL, _maxL);
  }

  /**
   * Saves the computed invariants into a binary file
   * @param _fName Output file name.
   */
  void SaveInvariants(const char *_fName) {
    std::ofstream outfile(_fName, std::ios_base::binary | std::ios_base::out);

    float temp;
    int dim = invariants_.size();
    outfile.write((char *)(&dim), sizeof(int));

    for (int i = 0; i < dim; ++i) {
      temp = invariants_[i];
      outfile.write((char *)(&temp), sizeof(float));
    }
  }

  /**
   * Public interface for accessing invariants.
   */
  T1D GetInvariants() { return invariants_; }

private:
  /**
   * Cuts elements from data object that fall outside unit sphere.
   */
  void NormalizeGrid() {
    T radius = (T)1 / scale_;
    T sqrRadius = radius * radius;

    T dn, dxy2, dxyz2;
    std::vector<T> dx2(dim_), dy2(dim_), dz2(dim_);

    for (int i = 0; i < dim_; ++i) {
      dn = (T)i - center[0];
      dx2[i] = dn * dn;
      dn = (T)i - center[1];
      dy2[i] = dn * dn;
      dn = (T)i - center[2];
      dz2[i] = dn * dn;
    }

    for (int x = 0; x < dim_; ++x) {
      for (int y = 0; y < dim_; ++y) {
        dxy2 = dx2[x] + dy2[y];
        for (int z = 0; z < dim_; ++z) {
          dxyz2 = dxy2 + dz2[z];
          if (dxyz2 > sqrRadius) {
            voxels_[(z * dim_ + y) * dim_ + x] = 0.0;
          }
        }
      }
    }
  }

  /**
   * Compute shift and scale to map input data to unit sphere.
   */
  void ComputeNormalization() {
    ScaledGeometricalMoments<T, T> gm(voxels_, dim_, dim_, dim_, 0.0, 0.0, 0.0,
                                      1.0);
    // 0'th order moments -> normalization, 1'st order center of gravity
    zeroMoment_ = gm.GetMoment(0, 0, 0);
    center[0] = gm.GetMoment(1, 0, 0) / zeroMoment_;
    center[1] = gm.GetMoment(0, 1, 0) / zeroMoment_;
    center[2] = gm.GetMoment(0, 0, 1) / zeroMoment_;

    T recScale = (T)2.0 * ComputeScale_RadiusVar();

    scale_ = (T)1 / recScale;
  }

  /**
   * Compute Zernike Moments.
   */
  void ComputeMoments() {
    gm_.Init(voxels_, dim_, dim_, dim_, center[0], center[1], center[2], scale_,
             order_);

    // Zernike moments
    zm_.Init(order_, gm_);
    zm_.Compute();
  }

  /**
   * Computes rotation invariant Moments as norm over Z_nl^m with m being the
   * running index.
   */
  void ComputeInvariants() {
    invariants_.clear();
    for (int n = 0; n < order_ + 1; ++n) {
      int l0 = n % 2, li = 0;
      for (int l = n % 2; l <= n; ++li, l += 2) {
        T sum = (T)0;

        for (int m = -l; m <= l; ++m) {
          ComplexT moment = zm_.GetMoment(n, l, m);
          sum += std::norm(moment);
        }
        invariants_.push_back(sqrt(sum));
      }
    }
  }

  /**
   * Computes the average distance between center and all larger-than zero
   * voxels.
   */
  T ComputeScale_RadiusVar() {
    T dn, dxy2, dxyz2;
    std::vector<T> dx2(dim_), dy2(dim_), dz2(dim_);

    for (int i = 0; i < dim_; ++i) {
      dn = (T)i - center[0];
      dx2[i] = dn * dn;
      dn = (T)i - center[1];
      dy2[i] = dn * dn;
      dn = (T)i - center[2];
      dz2[i] = dn * dn;
    }

    int d = dim_;
    T sum = 0.0;
    int nVoxels = 0;
    for (int x = 0; x < d; ++x) {
      for (int y = 0; y < d; ++y) {
        dxy2 = dx2[x] + dy2[y];
        for (int z = 0; z < d; ++z) {

          dxyz2 = dxy2 + dz2[z];
          if (std::fabs(voxels_[(z + d * y) * d + x]) > 0) {
            sum += dxyz2;
            nVoxels++;
          }
        }
      }
    }

    return (T)sqrt(sum / nVoxels);
  }

  /*
   * Reads cubic grid from a binary file.
   */
  T *ReadGrid(const char *_fname, int &dim_) {
    std::ifstream infile(_fname);
    if (!infile) {
      std::cerr << "Cannot open " << _fname << " for reading.\n";
      exit(-1);
    }

    std::vector<T> tempGrid;
    std::string line;
    TIn temp;

    // Read the grid values
    while (std::getline(infile, line)) {
      std::istringstream iss(line);
      if (iss >> temp) {
        tempGrid.push_back(static_cast<T>(temp));
      }
    }

    int d = tempGrid.size();
    double f = std::pow(static_cast<double>(d), 1.0 / 3.0);
    dim_ = static_cast<int>(std::floor(f + 0.5));

    T *result = new T[d];
    std::copy(tempGrid.begin(), tempGrid.end(), result);

    return result;
  }

private:
  int dim_;   // length of the edge of the voxel grid (which is a cube)
  int order_; // maximal order of the moments to be computed (max{n})

  T *voxels_;    // 1D array containing the voxels.
  T scale_;      // Scaling factor to map to unit shpere.
  T center[3];   // Center of gravity.
  T zeroMoment_; // Zero order moment.

  T1D invariants_; // Vector of invariants under SO(3)

  ZernikeMomentsT zm_;           // Actual Zernike Moments
  ScaledGeometricalMomentsT gm_; // Grid operations
};
