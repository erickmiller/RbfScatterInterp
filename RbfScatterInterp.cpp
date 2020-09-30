//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Library General Public
// License as published by the Free Software Foundation; either
// version 2 of the License, or (at your option) any later version.
// 
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Library General Public License for more details.
// 
// You should have received a copy of the GNU Library General Public
// License along with this library; if not, write to the
// Free Software Foundation, Inc., 59 Temple Place - Suite 330,
// Boston, MA  02111-1307, USA.
//
// Note: The simple 2d example class RbfScatterInterp.java by jplewis 
// was used as reference when writing this, and it is under this license 
// http://scribblethink.org/Work/PSD/RbfScatterInterp.java
// which dictates this license extends to any derivative works
//
#include <math.h>
#include <iostream>
#include <algorithm> 
#include <iostream> 
#include <vector> 
#define RBF_SCATTER_SQUARE(x) ((x)*(x)) 
using namespace std;

/**
*
* Author: Erick Miller
* Purpose: Radial Basis Function based interpolation in unlimited multi-dimensional vector space 
* Dependencies:  None. Standard C++ libraries. Algo and math are all self contained.
* Year written: 2004
*
* This is my C++ implementation of converting     
* Radial Basis Function scattered interpolation 
* I wrote this code entirely at home in evenings and weekends as a 
* fun chalenging exploratory exercise for the purpose creating a smooth 
* stable N dimensional interpolator where multi-dimensional vectors can be 
* driven and interpolated between using lower dimensional (like 1D) space
*
* Used a 2d version written in java as a reference implementation for this algorithm
* although we all know extending an algorithm like this from two dimensions
* into unlimited N dimensional space creates far more math and logic requirements
* Nonetheless, RbfScatterInterp.java, by jplewis was used as reference code
* and example code extensively when writing this library, thus deserves credit.  
* http://scribblethink.org/Work/PSD/
*
* This class is a standalone C++ library class that has
* no outside dependencies to non-standard external libraries, and 
* included is my ported/forked version of the matrix LU solve function. 
* Because it's standalone library, this can be easily included and used in 
* any program as an N-dimensional scattered data interpolation engine without
* the need to include multiple large dependencies as is normally required for 
* functionality such as this.
* 
*/

class RbfScatterInterp
{ 
    public:
        RbfScatterInterp( vector<vector<double> > &pts,  
            vector<double> &values, double sigma_mult, double reg);
        RbfScatterInterp(vector<vector<double> > &pts, 
                        vector<double> &values, double sigma_mult);
        RbfScatterInterp(vector<vector<double> > &pts,  
                                    vector<double> &values );
        virtual ~RbfScatterInterp();
        void 	setup(vector<vector<double> > &pts, 
                        vector<double> &values, double &width, double &reg);
        double 	interp(vector<double> &curr_pts);

    private: 
        double  distN(unsigned &a, unsigned &b, vector<vector<double> > &pts);
        double 	solve(vector<double> &M, vector<double> &b );
        double  computeSigma(vector<vector<double> > &pts ) ;
        unsigned				  _nd; 
        unsigned				  _npts;
        vector<vector<double> >	  _pts;
        vector<double>			  _w;
        double				  _msigma2;
};


RbfScatterInterp::RbfScatterInterp( vector<vector<double> > &pts,  
									vector<double> &values, 
									double sigma_mult, double reg) {
	double width = computeSigma(pts);
	width = width * sigma_mult ;
	setup(pts,values,width,reg);
}

RbfScatterInterp::RbfScatterInterp( vector<vector<double> > &pts,  
									vector<double> &values, 
									double sigma_mult) {
	double width = computeSigma(pts);
	width = width * sigma_mult ;
	double reg = 0.01 ;	
	setup(pts,values,width,reg);
}

RbfScatterInterp::RbfScatterInterp( vector<vector<double> > &pts,  
									vector<double> &values ) {
	double width = computeSigma(pts);
	double reg = 0.01 ;
	setup(pts,values,width,reg);
}
RbfScatterInterp::~RbfScatterInterp(){
	vector<vector<double> >	EMPTY_pts;
	vector<double>	  EMPTY_w;
	_pts = EMPTY_pts ;
	 _pts.clear();
	 _pts.resize(0);
	 _w = EMPTY_w;
	_w.clear();
	_w.resize(0);
};

/**
* TODO believe it is better if n-dimensional locations are 
* first mapped to near unity before solving this.
* 
* M(x,y) = exp(-(ax^2 + by^2) = exp(-ax^2)*exp(-by^2)
* [ M(0)  M(|p1 - p0|)  M(|p2 - p0|) ]   [w0x]     [val(0).x]
* [ M(|p0 - p1|)  M(1)  M(|p2 - p1|) ]   [w1x]  =  [val(1).x]
* [ M(|p0 - p2|)  M(|p1 - p2|)  M(0) ]   [w2x]     [val(2).x]
*
* Yhe above explanation is for a 2d case, and of course this same
* lienar algebra holds true for N-dimensional points and vectors
* --> N-D case simply splits each N into seperate 1D systems.
*/
void RbfScatterInterp::setup(vector<vector<double> > &pts,  
					vector<double> &values, double &width, double &reg){
		_npts = ((unsigned)pts.size());
		if(_npts == values.size() && _npts>0 ){ 
		_pts = pts;
		_nd = ((unsigned)pts[0].size() );
		////
		// TODO: Make the width an array of doubles
		//   accessable from within the inner loop:
		//
		width = 1.0 / width;
		_msigma2 = -(width * width);
		vector<double> M ;
		M.resize( RBF_SCATTER_SQUARE(_npts) );
		unsigned mm=0;
		double nd=0.0;	
		for( unsigned r=0; r < _npts; r++ ){
			for( unsigned c=0; c < _npts; c++, mm++ ){
				nd = distN(r, c, pts);
				M[mm] = exp(_msigma2 * nd);
				if(r==c)
				{
					/////
					// The following ridge-regression / matrix regularization at the 
					// diagonal introduces tons of good stability to the interpolant 
                    // when space is high dimensional, overlapping or not well defined,  
					// with the very small cost of just the slightest delta inaccuracy:
					//
					// TODO: Add this as a per constraint pt array argument:
					//
					M[mm] = M[mm] + reg ;
				}
			}
		}
		////
		// TODO: Warn about this determinant:
		//
		double determinant = solve(M, values);
		_w = values;
	}
} //setup

double RbfScatterInterp::interp(vector<double> &curr_pts){
	double sum=0.0, sum_squared=0.0;
	if( _nd == curr_pts.size() && _nd>0 ){
		for( unsigned i=0; i < _npts; i++ ){
			sum_squared=0.0;
			for( unsigned n=0;n<curr_pts.size();n++ ){
				sum_squared += RBF_SCATTER_SQUARE( (curr_pts[n] - _pts[i][n]) );
			}
			sum += (_w[i] * exp(_msigma2 * sum_squared));//expensive
		}
	}
	return sum;
} //interp

/**
* Return distance (w/ out square root) between n-dimensional points a,b
*/
double RbfScatterInterp::distN(unsigned &a, unsigned &b, 
								  vector<vector<double> > &pts ){
	double sum_squared = 0.0;
	for( unsigned n=0; n < pts[a].size(); n++ ){
		sum_squared += RBF_SCATTER_SQUARE( pts[a][n] - pts[b][n] );
	}
	return( sum_squared );
}

/**
* 	Solve an nxn system of linear equations Mx=b
* 	using Gaussian elimination with partial pivoting.
* 	leave solution x in b array (destroying original M and 
* 	b in the process) Returns determinant. 
*/
double RbfScatterInterp::solve(vector<double> &M, vector<double> &b)
{
#  define swap(a,b,t) {t=a; a=b; b=t;}
#  define a(i,j) M[(i)*n+(j)]

   int i,j,k, n = ((int)b.size());
   double max,t,det,sum,pivot;	/* keep these double */
   /*---------- forward elimination ----------*/
   det = 1.0;
   for (i=0; i<n; i++) {		/* eliminate in column i */
      max = -1.;
      for (k=i; k<n; k++)		/* find pivot for column i */
         if (fabs(a(k,i))>max) {
            max = fabs(a(k,i));
            j = k;
         }
      if (max<=0.0) return(0.0);/* if no nonzero pivot, PUNT */
      if (j!=i) {				/* swap rows i and j */
         for (k=i; k<n; k++)
            swap(a(i,k),a(j,k),t);
         det = -det;
         swap(b[i],b[j],t);		/* swap elements of column vector */
      }
      pivot = a(i,i);
      det *= pivot;
      for (k=i+1; k<n; k++)		/* only do elems to right of pivot */
         a(i,k) /= pivot;
     /* we know that a(i,i) will be set to 1, so why bother to do it? */
      b[i] /= pivot;
      for (j=i+1; j<n; j++) {	/* eliminate in rows below i */
         t = a(j,i);			/* we're gonna zero this guy */
         for (k=i+1; k<n; k++)	/* subtract scaled row i from row j */
            a(j,k) -= a(i,k)*t;	/* (ignore k<=i, we know they're 0) */
         b[j] -= b[i]*t;
      }
   }
   /*---------- back substitution ----------*/
   for (i=n-1; i>=0; i--) {		/* solve for x[i] (put it in b[i]) */
      sum = b[i];
      for (k=i+1; k<n; k++)		/* really a(i,k)*x[k] */
         sum -= a(i,k)*b[k];
      b[i] = sum;
   }
   return(det);

#  undef swap
#  undef a
} /*solve*/

/**
* Uses n-dimensional vector lengths to derive a
* default one-dimensional width from n-space, using the 
* difference between the highest and lowest n-dimensional 
* lengths as the width for sigma in the psd equation.
*/
double RbfScatterInterp::computeSigma(vector<vector<double> > &pts)
{
	unsigned num = pts.size();
	double sum_squared = 0.0;
	vector<double> his(num);
	vector<double> los(num);
	double curr,hi,lo;
	unsigned i=0;
	for(i=0;i<num;i++){
		curr=0.0; hi=-9e99; lo=9e99;
		for(unsigned p=0;p<num;p++){
			sum_squared = 0.0;
			unsigned nd = pts[p].size();
			for(unsigned n=0;n<nd;n++){
				sum_squared += RBF_SCATTER_SQUARE( pts[i][n]-pts[p][n] );
			} 
			curr = sqrt( sum_squared ); 
			if(curr>hi){ hi=curr; his[i]=hi; }
			if(curr<lo){ lo=curr; los[i]=lo; }
		}
	}
	hi=-9e99; lo=9e99;
	for( i=0; i<num; i++ ){
		if(his[i]>hi)  
			hi = his[i];
		if(los[i]<lo)  
			lo = los[i];
	}
	return ( hi-lo );
}

#undef RBF_SCATTER_SQUARE