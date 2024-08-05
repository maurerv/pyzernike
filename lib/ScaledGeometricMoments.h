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

#include <vector>

/**
 * Class for computing the scaled, pre-integrated geometrical moments.
 */
template <class VoxelT, class MomentT> class ScaledGeometricalMoments {
public:
  typedef MomentT T;                      // the moment type
  typedef std::vector<T> T1D;             // vector scalar type
  typedef std::vector<T1D> T2D;           // 2D array scalar type
  typedef std::vector<T2D> T3D;           // 3D array scalar type
  typedef std::vector<double> Double1D;   // vector scalar type
  typedef std::vector<Double1D> Double2D; // vector scalar type
  typedef typename T1D::iterator T1DIter;

  ScaledGeometricalMoments() {}

  /**
   * Constructor.
   *
   * @param _voxels Input voxel grid.
   * @param _xDim x-dimension of the input voxel grid.
   * @param _yDim z-dimension of the input voxel grid.
   * @param _yDim y-dimension of the input voxel grid.
   * @param _xCOG x-coord of the center of gravity.
   * @param _yCOG y-coord of the center of gravity.
   * @param _zCOG z-coord of the center of gravity.
   * @param _scale Scaling factor.
   * @param _maxOrderMaximal order to compute moments for.
   */
  ScaledGeometricalMoments(
      const VoxelT *_voxels, /**< input voxel grid */
      int _xDim,             /**< x-dimension of the input voxel grid */
      int _yDim,             /**< y-dimension of the input voxel grid */
      int _zDim,             /**< z-dimension of the input voxel grid */
      double _xCOG,          /**< x-coord of the center of gravity */
      double _yCOG,          /**< y-coord of the center of gravity */
      double _zCOG,          /**< z-coord of the center of gravity */
      double _scale,         /**< scaling factor */
      int _maxOrder = 1      /**< maximal order to compute moments for */
  ) {
    Init(_voxels, _xDim, _yDim, _zDim, _xCOG, _yCOG, _zCOG, _scale, _maxOrder);
  }

  /**
   * The init function used by the contructors.
   */
  void Init(const VoxelT *_voxels, int _xDim, int _yDim, int _zDim,
            double _xCOG, double _yCOG, double _zCOG, double _scale,
            int _maxOrder = 1) {
    xDim_ = _xDim;
    yDim_ = _yDim;
    zDim_ = _zDim;

    maxOrder_ = _maxOrder;

    size_t totalSize = xDim_ * yDim_ * zDim_;
    voxels_.resize(totalSize);
    for (int i = 0; i < totalSize; ++i) {
      voxels_[i] = _voxels[i];
    }

    moments_.resize(maxOrder_ + 1);
    for (int i = 0; i <= maxOrder_; ++i) {
      moments_[i].resize(maxOrder_ - i + 1);
      for (int j = 0; j <= maxOrder_ - i; ++j) {
        moments_[i][j].resize(maxOrder_ - i - j + 1);
      }
    }

    ComputeSamples(_xCOG, _yCOG, _zCOG, _scale);

    Compute();
  }

  /**
   *  Access Geometrical Moments.
   *
   * @param _n Order along x.
   * @param _l Order along y.
   * @param _m Order along z.
   * @return Complex Moment.
   */
  T GetMoment(int _i, int _j, int _k) { return moments_[_i][_j][_k]; }

private:
  void Compute();

  void ComputeSamples(double _xCOG, double _yCOG, double _zCOG, double _scale) {
    samples_.resize(3); // 3 dimensions

    int dim[3] = {xDim_, yDim_, zDim_};
    double min[3] = {(-_xCOG) * _scale, (-_yCOG) * _scale, (-_zCOG) * _scale};
    for (int i = 0; i < 3; ++i) {
      samples_[i].resize(dim[i] + 1);
      for (int j = 0; j <= dim[i]; ++j) {
        samples_[i][j] = min[i] + j * _scale;
      }
    }
  }

  void ComputeDiffFunction(T1DIter _iter, T1DIter _diffIter, int _dim) {
    _diffIter[0] = -_iter[0];
    for (int i = 1; i < _dim; ++i) {
      _diffIter[i] = _iter[i - 1] - _iter[i];
    }
    _diffIter[_dim] = _iter[_dim - 1];
  }

  T Multiply(T1DIter _diffIter, T1DIter _sampleIter, int _dim) {
    T sum = 0;
    for (int i = 0; i < _dim; ++i) {
      _diffIter[i] *= _sampleIter[i];
      sum += _diffIter[i];
    }
    return sum;
  }

private:
  int maxOrder_;           // maximal order of the moments
  int xDim_, yDim_, zDim_; // Data dimensions

  T1D voxels_;  // array containing the voxel grid
  T2D samples_; // samples of the scaled and translated grid in x, y, z
  T3D moments_; // array containing the cumulative moments
};

template <class VoxelT, class MomentT>
void ScaledGeometricalMoments<VoxelT, MomentT>::Compute() {
  int arrayDim = zDim_;
  int layerDim = yDim_ * zDim_;

  int diffArrayDim = zDim_ + 1;
  int diffLayerDim = (yDim_ + 1) * zDim_;
  int diffGridDim = (xDim_ + 1) * layerDim;

  T1D diffGrid(diffGridDim);
  T1D diffLayer(diffLayerDim);
  T1D diffArray(diffArrayDim);

  T1D layer(layerDim);
  T1D array(arrayDim);
  T moment;

  typename T1D::iterator iter = voxels_.begin();
  typename T1D::iterator diffIter = diffGrid.begin();

  // generate the diff version of the voxel grid in x direction
  for (int x = 0; x < layerDim; ++x) {
    ComputeDiffFunction(iter, diffIter, xDim_);

    iter += xDim_;
    diffIter += xDim_ + 1;
  }

  for (int i = 0; i <= maxOrder_; ++i) {
    diffIter = diffGrid.begin();
    for (int p = 0; p < layerDim; ++p) {
      // multiply the diff function with the sample values
      T1DIter sampleIter(samples_[0].begin());
      layer[p] = Multiply(diffIter, sampleIter, xDim_ + 1);

      diffIter += xDim_ + 1;
    }

    iter = layer.begin();
    diffIter = diffLayer.begin();
    for (int y = 0; y < arrayDim; ++y) {
      ComputeDiffFunction(iter, diffIter, yDim_);

      iter += yDim_;
      diffIter += yDim_ + 1;
    }

    for (int j = 0; j < maxOrder_ + 1 - i; ++j) {
      diffIter = diffLayer.begin();
      for (int p = 0; p < arrayDim; ++p) {
        T1DIter sampleIter(samples_[1].begin());
        array[p] = Multiply(diffIter, sampleIter, yDim_ + 1);

        diffIter += yDim_ + 1;
      }

      iter = array.begin();
      diffIter = diffArray.begin();
      ComputeDiffFunction(iter, diffIter, zDim_);

      for (int k = 0; k < maxOrder_ + 1 - i - j; ++k) {
        T1DIter sampleIter(samples_[2].begin());

        moment = Multiply(diffIter, sampleIter, zDim_ + 1);
        moments_[i][j][k] = moment / ((1 + i) * (1 + j) * (1 + k));
      }
    }
  }
}
