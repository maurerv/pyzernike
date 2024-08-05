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
    conference = {The 8th ACM Symposium on Solid Modeling and Applications, June 16-20, Seattle, WA}
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

#ifndef ZERNIKEMOMENTS_H
#define ZERNIKEMOMENTS_H

#pragma warning (disable : 4267)

// ----- local program includes -----
#include "Factorial.h"
#include "Binomial.h"
#include "ScaledGeometricMoments.h"

// ----- std includes -----
#include <complex>
#include <cstdlib>
#include <set>
#include <ios>

#define PI 3.141592653589793

/**
 * Struct representing a complex coefficient of a moment
 * of order (p_,q_,r_)
 */
template<class T>
struct ComplexCoeff
{
    typedef     std::complex<T>     ValueT;

    ComplexCoeff (int, int, int, const ValueT&);
    ComplexCoeff (const ComplexCoeff<T>& _cc);
    ComplexCoeff ();

    int                 p_, q_, r_;
    ValueT     value_;
};

/**
 * Class representing the Zernike moments
 */
template<class VoxelT, class MomentT>
class ZernikeMoments
{
public:
    // ---- public typedefs ----
    typedef MomentT             T;
    typedef vector<T>           T1D;        // vector of scalar type
    typedef vector<T1D>         T2D;        // 2D array of scalar type
    typedef vector<T2D>         T3D;        // 3D array of scalar type
    typedef vector<T3D>         T4D;        // 3D array of scalar type

    typedef std::complex<T>                      ComplexT;       // complex type
    typedef vector<vector<vector<ComplexT> > >   ComplexT3D;     // 3D array of complex type

    typedef ComplexCoeff<T>                      ComplexCoeffT;
    typedef vector<vector<vector<vector<ComplexCoeffT> > > >    ComplexCoeffT4D;


public:
    // ---- public member functions ----
    ZernikeMoments (int _order, ScaledGeometricalMoments<VoxelT,MomentT>& _gm);
    ZernikeMoments ();

    /* Todo: Init() doesn't depend on the geomentrical moments.
     * Removing the  parameter and moving it to Compute() will make the class more efficient to reuse to
     * compute ZM for several objects at once.
     */
    void Init (int _order, ScaledGeometricalMoments<VoxelT,MomentT>& _gm);
    void Compute ();

    ComplexT GetMoment (int _n, int _l, int _m);

    // ---- debug functions/arguments ----
    void Reconstruct (ComplexT3D&   _grid,                // grid containing the reconstructed function
                      T             _xCOG,                // center of gravity
                      T             _yCOG,
                      T             _zCOG,
                      T             _scale,               // scaling factor to map into unit ball
                      int           _minN = 0,            // min value for n freq index
                      int           _maxN = 100,          // min value for n freq index
                      int           _minL = 0,            // min value for l freq index
                      int           _maxL = 100);         // max value for l freq index

    void NormalizeGridValues (ComplexT3D& _grid);
    void CheckOrthonormality (int _n1, int _l1, int _m1, int _n2, int _l2, int _m2);

private:
    // ---- private member functions ----
    void ComputeCs ();
    void ComputeQs ();
    void ComputeGCoefficients ();

    // ---- private attributes -----
    ComplexCoeffT4D     gCoeffs_;           // coefficients of the geometric moments
    ComplexT3D          zernikeMoments_;    // nomen est omen
    T3D                 qs_;                // q coefficients (radial polynomial normalization)
    T2D                 cs_;                // c coefficients (harmonic polynomial normalization)

    ScaledGeometricalMoments<VoxelT,MomentT> gm_;
    int                 order_;             // := max{n} according to indexing of Zernike polynomials

    // ---- debug functions/arguments ----
    void PrintGrid (ComplexT3D& _grid);
    T EvalMonomialIntegral (int _p, int _q, int _r, int _dim);
};


// Moved from ZernikeMoments.inl
template<class VoxelT, class MomentT>
inline typename ZernikeMoments<VoxelT,MomentT>::ComplexT
ZernikeMoments<VoxelT,MomentT>::GetMoment (int _n, int _l, int _m)
{
    if (_m >= 0)
    {
        return zernikeMoments_[_n][_l/2][_m];
    }
    else
    {
        T sign;
        if (_m%2)
        {
            sign = (T)(-1);
        }
        else
        {
            sign = (T)1;
        }
        return sign * std::conj (zernikeMoments_[_n][_l/2][abs(_m)]);
    }
}


// Moved from ZernikeMoments.cpp
/**
 * Copy constructor
 */
template<class T>
ComplexCoeff<T>::ComplexCoeff (const ComplexCoeff<T>& _cc) :
    p_ (_cc.p_), q_ (_cc.q_), r_ (_cc.r_), value_ (_cc.value_)
{
}

/**
 * Constructor with scalar args
 */
template<class T>
ComplexCoeff<T>::ComplexCoeff (int _p, int _q, int _r, const ValueT& _value) :
    p_ (_p), q_ (_q), r_ (_r), value_ (_value)
{
}

/**
 * Default constructor
 */
template<class T>
ComplexCoeff<T>::ComplexCoeff () :
    p_ (0), q_ (0), r_ (0)
{
}


// ---------- implementation of ZernikeMoments class -------------
template<class VoxelT, class MomentT>
ZernikeMoments<VoxelT,MomentT>::ZernikeMoments (int _order, ScaledGeometricalMoments<VoxelT,MomentT>& _gm)
{
    Init (_order, _gm);
}

template<class VoxelT, class MomentT>
ZernikeMoments<VoxelT,MomentT>::ZernikeMoments () :
    order_ (0)
{
}

/**
 * Computes all coefficients that are input data independent
 */
template<class VoxelT, class MomentT>
void ZernikeMoments<VoxelT,MomentT>::Init (int _order, ScaledGeometricalMoments<VoxelT,MomentT>& _gm)
{
    /* Todo: Init() doesn't depend on the geomentrical moments.
     * Removing the  parameter and moving it to Compute() will make the class more efficient to reuse to
     * compute ZM for several objects at once.
     */
    gm_ = _gm;
    order_ = _order;

    ComputeCs ();
    ComputeQs ();
    ComputeGCoefficients ();
}


/**
 * Computes all the normalizing factors $c_l^m$ for harmonic polynomials e
 */
template<class VoxelT, class MomentT>
void ZernikeMoments<VoxelT,MomentT>::ComputeCs ()
{
    /*
     indexing:
       l goes from 0 to n
       m goes from -l to l, in fact from 0 to l, since c(l,-m) = c (l,m)
    */

    cs_.resize (order_ + 1);

    for (int l=0; l<=order_; ++l)
    {
        cs_[l].resize (l + 1);
        for (int m=0; m<=l; ++m)
        {
            T n_sqrt = ((T)2 * l + (T)1) *
                       Factorial<T>::Get (l + 1, l + m);
            T d_sqrt = Factorial<T>::Get (l - m + 1, l);

            cs_[l][m] = sqrt (n_sqrt / d_sqrt);
        }
    }
}

/**
 * Computes all coefficients q for orthonormalization of radial polynomials
 * in Zernike polynomials.
 */
template<class VoxelT, class MomentT>
void ZernikeMoments<VoxelT,MomentT>::ComputeQs ()
{
    /*
     indexing:
       n goes 0..order_
       l goes 0..n, so that n-l is even
       mu goes 0..(n-l)/2
    */

    qs_.resize (order_ + 1);            // there is order_ + 1 n's

    for (int n=0; n<=order_; ++n)
    {
        qs_[n].resize (n / 2 + 1);      // there is floor(n/2) + 1 l's

        int l0 = n % 2;
        for (int l=l0; l<=n; l+=2)
        {
            int k = (n-l)/2;

            qs_[n][l/2].resize (k + 1);   // there is k+1 mu's

            for (int mu=0; mu<=k; ++mu)
            {
                T nom = Binomial<T>::Get (2*k, k) *     // nominator of straight part
                        Binomial<T>::Get (k, mu) *
                        Binomial<T>::Get (2 * (k + l + mu) + 1, 2 * k);

                if ((k+mu) % 2)
                {
                    nom *= (T)(-1);
                }

                T den = std::pow ((T)2, (T)(2*k)) *     // denominator of straight part
                        Binomial<T>::Get (k + l + mu, k);

                T n_sqrt = (T)(2 * l + 4 * k + 3);      // nominator of sqrt part
                T d_sqrt = (T)3;                        // denominator of sqrt part

                qs_[n][l/2][mu] =  nom / den * sqrt (n_sqrt / d_sqrt);
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
template<class VoxelT, class MomentT>
void ZernikeMoments<VoxelT,MomentT>::ComputeGCoefficients ()
{
//DD
    int countCoeffs = 0;
//DD
    gCoeffs_.resize (order_+1);
    for (int n=0; n<=order_; ++n)
    {
        gCoeffs_[n].resize (n/2+1);
        int li = 0, l0 = n%2;
        for (int l=l0; l<=n; ++li, l+=2)
        {
            gCoeffs_[n][li].resize (l+1);
            for (int m=0; m<=l; ++m)
            {
                T w = cs_[l][m] / std::pow ((T)2, (T)m);

                int k= (n-l)/2;
                for (int nu=0; nu<=k; ++nu)
                {
                    T w_Nu = w * qs_[n][li][nu];
                    for (int alpha=0; alpha<=nu; ++alpha)
                    {
                        T w_NuA = w_Nu * Binomial<T>::Get (nu, alpha);
                        for (int beta=0; beta<=nu-alpha; ++beta)
                        {
                            T w_NuAB = w_NuA * Binomial<T>::Get (nu-alpha, beta);
                            for (int p=0; p<=m; ++p)
                            {
                                T w_NuABP = w_NuAB * Binomial<T>::Get (m, p);
                                for (int mu=0; mu<=(l-m)/2; ++mu)
                                {
                                    T w_NuABPMu = w_NuABP *
                                                  Binomial<T>::Get (l, mu) *
                                                  Binomial<T>::Get (l-mu, m+mu) /
                                                  (T)pow (2.0, (double)(2*mu));
                                    for (int q=0; q<=mu; ++q)
                                    {
                                        // the absolute value of the coefficient
                                        T w_NuABPMuQ = w_NuABPMu * Binomial<T>::Get (mu, q);

                                        // the sign
                                        if ((m-p+mu)%2)
                                        {
                                            w_NuABPMuQ *= T(-1);
                                        }

                                        // * i^p
                                        int rest = p % 4;
                                        ComplexT c;
                                        switch (rest)
                                        {
                                        case 0: c = ComplexT (w_NuABPMuQ, (T)0); break;
                                        case 1: c = ComplexT ((T)0, w_NuABPMuQ); break;
                                        case 2: c = ComplexT ((T)(-1) * w_NuABPMuQ, (T)0); break;
                                        case 3: c = ComplexT ((T)0, (T)(-1) * w_NuABPMuQ); break;
                                        }

                                        // determination of the order of according moment
                                        int z_i = l-m+2*(nu-alpha-beta-mu);
                                        int y_i = 2*(mu-q+beta)+m-p;
                                        int x_i = 2*q+p+2*alpha;
//DD
                                        //std::cout << x_i << " " << y_i << " " << z_i;
                                        //std::cout << "\t" << n << " " << l << " " << m;
                                        //std::cout << "\t" << c.real () << " " << c.imag () << std::endl;
//DD
                                        ComplexCoeffT cc (x_i, y_i, z_i, c);
                                        gCoeffs_[n][li][m].push_back (cc);
//DD
                                        countCoeffs++;
//DD
                                    } // q
                                } // mu
                            } // p
                        } // beta
                    } // alpha
                } // nu
            } // m
        } // l
    } // n
//DD
    std::cout << countCoeffs << std::endl;
//DD
}

/**
 * Computes the Zernike moments. This computation is data dependent
 * and has to be performed for each new object and/or transformation.
 */
template<class VoxelT, class MomentT>
void ZernikeMoments<VoxelT,MomentT>::Compute ()
{
    // geometrical moments have to be computed first
    if (!order_)
    {
        std::cerr << "ZernikeMoments<VoxelT,MomentT>::ComputeZernikeMoments (): attempting to \
                     compute Zernike moments without setting valid geometrical \
                     moments first. \nExiting...\n";
        exit (-1);
    }

    /*
     indexing:
       n goes 0..order_
       l goes 0..n so that n-l is even
       m goes -l..l
    */

    T nullMoment;

    zernikeMoments_.resize (order_ + 1);
    for (int n=0; n <= order_; ++n)
    {
        zernikeMoments_[n].resize (n/2 + 1);

        int l0 = n % 2, li = 0;
        for (int l = l0; l<=n; ++li, l+=2)
        {
            zernikeMoments_[n][li].resize (l + 1);
            for (int m=0; m<=l; ++m)
            {
                // Zernike moment of according indices [nlm]
                ComplexT zm ((T)0, (T)0);

                int nCoeffs = (int)gCoeffs_[n][li][m].size ();
                for (int i=0; i<nCoeffs; ++i)
                {
                    ComplexCoeffT cc = gCoeffs_[n][li][m][i];
                    //T scale = gm_.GetScale ();
                    //T fact = std::pow (scale, cc.p_+cc.q_+cc.r_+3);

                    //zm +=  std::conj (cc.value_) * gm_.GetMoment(cc.p_, cc.q_, cc.r_) * fact;
                    zm +=  std::conj (cc.value_) * gm_.GetMoment(cc.p_, cc.q_, cc.r_);
                }

                zm *= (T)(3.0 / (4.0 * PI));
                if (n ==0 && l == 0 && m == 0)
                {
                    // FixMe Unused variable! What is it for?
                    nullMoment = zm.real ();
                }
                zernikeMoments_[n][li][m] = zm;

                //std::cout << "zernike moment[nlm]: " << n << "\t" << l << "\t" << m << "\t" << zernikeMoments_[n][li][m] << "\n";
            }
        }
    }
}


/**
 * The function previously encoded as complex valued Zernike
 * moments, is reconstructed. _grid is the output grid containing
 * the reconstructed function.
 */
template<class VoxelT, class MomentT>
void ZernikeMoments<VoxelT,MomentT>::Reconstruct (ComplexT3D& _grid, T _xCOG, T _yCOG, T _zCOG, T _scale, int _minN, int _maxN, int _minL, int _maxL)
{
    int dimX = _grid.size ();
    int dimY = _grid[0].size ();
    int dimZ = _grid[0][0].size ();

    //translation
    T vx = _xCOG;
    T vy = _yCOG;
    T vz = _zCOG;

    T dx, dy, dz;
    T dx2, dxy2, dxyz2;

    if (_maxN == 100) _maxN = order_;

    // origin is at the grid center, all voxels are projected onto the unit sphere
    for (int x=0; x<dimX; ++x) {
        std::cout << x << ". layer being processed\n";
        dx = ((T)x-vx) * _scale;
        dx2 = dx * dx;
        for (int y=0; y<dimY; ++y) {
            dy = ((T)y-vy) * _scale;
            dxy2 = dx2 + dy * dy;
            if (dxyz2 > 1.0) continue;

            for (int z=0; z<dimZ; ++z) {
                T dz = ((T)z - vz) * _scale;
                dxyz2 = dxy2 + dz * dz;
                if (dxyz2 > 1.0) continue;

                // function value at point
                ComplexT fVal = (0, 0);
                for (int n=_minN; n<=_maxN; ++n) {
                    int maxK = n/2;
                    for (int k=0; k<=maxK; ++k) {
                        for (int nu=0; nu<=k; ++nu) {
                            // check whether l is within bounds
                            int l=n-2*k;
                            if (l < _minL || l > _maxL) continue;

                            for (int m=-l; m<=l; ++m) {
                                // zernike polynomial evaluated at point
                                ComplexT zp (0, 0);

                                int absM = std::abs (m);

                                int nCoeffs = gCoeffs_[n][l/2][absM].size ();
                                for (int i=0; i<nCoeffs; ++i) {
                                    ComplexCoeffT cc = gCoeffs_[n][l/2][absM][i];
                                    ComplexT cvalue = cc.value_;

                                    // conjugate if m negative and take care of sign
                                    if (m<0) {
                                        cvalue = std::conj (cvalue);
                                        if (m % 2) cvalue *= (T)(-1);
                                    }

                                    zp += cvalue *
                                          std::pow (dx,(T)cc.p_) *
                                          std::pow (dy,(T)cc.q_) *
                                          std::pow (dz,(T)cc.r_);
                                }

                                fVal += zp * GetMoment(n, l, m);
                            }
                        }
                    }
                }
                _grid[x][y][z] = fVal;
            }
        }
    }

    //NormalizeGridValues (_grid);
}


template<class VoxelT, class MomentT>
void ZernikeMoments<VoxelT,MomentT>::NormalizeGridValues (ComplexT3D& _grid)
{
    int xD = _grid.size ();
    int yD = _grid[0].size ();
    int zD = _grid[0][0].size ();

    T max = (T)0;
    for (int k=0; k<zD; ++k)
    {
        for (int j=0; j<yD; ++j)
        {
            for (int i=0; i<xD; ++i)
            {
                if (_grid[i][j][k].real () > max)
                {
                    max = _grid[i][j][k].real ();
                }
            }
        }
    }

    std::cout << "\nMaximal value in grid: " << max << "\n";

    T invMax = (T)1 / max;

    for (int k=0; k<zD; ++k)
    {
        for (int j=0; j<yD; ++j)
        {
            for (int i=0; i<xD; ++i)
            {
                _grid[i][j][k] *= invMax;
            }
        }
    }
}


template<class VoxelT, class MomentT>
void ZernikeMoments<VoxelT,MomentT>::PrintGrid (ComplexT3D& _grid)
{
    int xD = _grid.size ();
    int yD = _grid[0].size ();
    int zD = _grid[0][0].size ();

    std::cout.setf (std::ios_base::scientific, std::ios_base::floatfield);

    T max = (T)0;
    for (int k=0; k<zD; ++k)
    {
        for (int j=0; j<yD; ++j)
        {
            for (int i=0; i<xD; ++i)
            {
                if (fabs (_grid[i][j][k].real ()) > max)
                {
                    max = fabs (_grid[i][j][k].real ());
                }
            }
        }
    }

    for (int k=0; k<zD; ++k)
    {
        std::cout << k << ". layer:\n";
        for (int j=0; j<yD; ++j)
        {
            for (int i=0; i<xD; ++i)
            {
                //std::cout << _grid[i][j][k].real () / max << "\t";
                std::cout << _grid[i][j][k].real () << "\t";
            }
            std::cout << "\n";
        }
        std::cout << "\n";
    }
    std::cout.setf (std::ios_base::fmtflags (0), std::ios_base::floatfield);
}

template<class VoxelT, class MomentT>
void ZernikeMoments<VoxelT,MomentT>::CheckOrthonormality (int _n1, int _l1, int _m1, int _n2, int _l2, int _m2)
{
    int li1 = _l1/2;
    int li2 = _l2/2;
    int dim = 64;

    // the total sum of the scalar product
    ComplexT sum ((T)0, (T)0);

    int nCoeffs1 = (int)gCoeffs_[_n1][li1][_m1].size ();
    int nCoeffs2 = (int)gCoeffs_[_n2][li2][_m2].size ();

    for (int i=0; i<nCoeffs1; ++i)
    {
        ComplexCoeffT cc1 = gCoeffs_[_n1][li1][_m1][i];
        for (int j=0; j<nCoeffs2; ++j)
        {
            ComplexCoeffT cc2 = gCoeffs_[_n2][li2][_m2][j];

            T temp = (T)0;

            int p = cc1.p_ + cc2.p_;
            int q = cc1.q_ + cc2.q_;
            int r = cc1.r_ + cc2.r_;

            sum +=  cc1.value_ *
                    std::conj (cc2.value_) *
                    EvalMonomialIntegral (p, q, r, dim);
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
template<class VoxelT, class MomentT>
MomentT ZernikeMoments<VoxelT,MomentT>::EvalMonomialIntegral (int _p, int _q, int _r, int _dim)
{
    T radius = (T)(_dim-1)/(T)2;
    T scale =  std::pow ((T)1/radius, 3);
    T center = (T)(_dim-1)/(T)2;

    T result = (T)0;
    T point[3];

    for (int x=0; x<_dim; ++x)
    {
        point[0] = ((T)x-center) / radius;
        for (int y=0; y<_dim; ++y)
        {
            point[1] = ((T)y-center) / radius;
            for (int z=0; z<_dim; ++z)
            {
                point[2] = ((T)z-center) / radius;

                if (point[0]*point[0] + point[1]*point[1] + point[2]*point[2] > (T)1)
                {
                    continue;
                }

                result += std::pow (point[0],(T)_p) *
                          std::pow (point[1],(T)_q) *
                          std::pow (point[2],(T)_r);
            }
        }
    }

    result *=  (T)(3.0 / (4.0 * PI)) * scale;
    return result;
}


#endif
