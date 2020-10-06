// RbfScatterInterp.java jplewis - psd scatterinterp engine
//
// * TODO believe it is better if x,y locations are first mapped
// * to near unity before solving this.
//
// modified
// febst scruff	

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
// contact info:  zilla@computer.org

package ZS.Grid;

import VisualNumerics.math.DoubleMatrix;

import zlib.*;

// uses double because we are using jnl.*
final public class RbfScatterInterp implements ScatterInterpSparse2
{
  int		_verbose = 0;
  int		_npts;
  double[][]	_pts;
  double[]	_w;
  double	_msigma2;

  /**
   * TODO also do 'wendland' function:
   * gnuplot> plot [x=0:1] (((1-x)**4)*(4*x+1)) 
   */
  public RbfScatterInterp(double[][] pts, double[] values, double width)
  {
    setup(pts, values, width);
  } //constructor


  /**
   */
  public RbfScatterInterp(float[][] fpts, float[] fvalues, float width)
  {
    double[][] pts = zlib.toDouble(fpts);
    double[] values = zlib.toDouble(fvalues);
    setup(pts, values, (double)width);
  } //constructor

  /**
   * TODO believe it is better if x,y locations are first mapped
   * to near unity before solving this.
   * 
   * M(x,y) = exp(-(ax^2 + by^2) = exp(-ax^2)*exp(-by^2)
   * [ M(0)  M(|p1 - p0|)  M(|p2 - p0|) ]   [w0x]     [val(0).x]
   * [ M(|p0 - p1|)  M(1)  M(|p2 - p1|) ]   [w1x]  =  [val(1).x]
   * [ M(|p0 - p2|)  M(|p1 - p2|)  M(0) ]   [w2x]     [val(2).x]
   *
   * 2D case splits into 2 1D systems, one for x-displacement, one for y.
   */
  private void setup(double[][] pts, double[] values, double width)
  {
    _npts = pts.length;
    zliberror.assert(_npts == values.length);
    _pts = zlib.cloneArray(pts);
    //_values = (double[])values.clone();
    width = 1.f / width;
    _msigma2 = -(width * width);

    double[][] M = new double[_npts][_npts];
    for( int r=0; r < _npts; r++ ) {
      for( int c=0; c < _npts; c++ ) {
	double d2 = dist2(r, c, pts);
	M[r][c] = Math.exp(_msigma2 * d2);
      }
    }

    try {
      // TODO: warn about determinant
      _w = DoubleMatrix.solve(M, values);
    }
    catch(Exception x) { zliberror.die(x); }

    if (_verbose > 0) {
      for( int i=0; i < _npts; i++ ) {
	double v = interp((float)pts[i][0], (float)pts[i][1]);
	System.out.println("wanted " + values[i] +", got " + v);
      }
    } 
  } //setup

  //----------------------------------------------------------------

  /**
   * rbf = sum w_k*G(| pt - pt[k] |)
   * msigma2 is -(_sigma*_sigma) 
   */
  public final float interp(float xf, float yf)
  {
    double sum = 0.f;
    double x = xf;
    double y = yf;
    int npts = _npts;

    for( int i=0; i < npts; i++ ) {
      double dx = x - _pts[i][0];
      double dy = y - _pts[i][1];
      double dist2 = dx*dx + dy*dy;
      sum += (_w[i] * Math.exp(_msigma2 * dist2));	// expensive
    }

    return (float)sum;
  } //interp


  //----------------------------------------------------------------


  /** return distance between points a,b
   */
  private static final double dist2(int a, int b, double[][] pts)
  {
    double dx = pts[a][0] - pts[b][0];
    double dy = pts[a][1] - pts[b][1];
    return dx*dx + dy*dy;
  }


  //----------------------------------------------------------------

  // test
  public static void main(String[] cmdline)
  {
    // 0          1
    //      .5
    // 1          ?

    float[][] pts = new float[][]
    {{0.f,0.f}, {1.f,0.f}, {0.f,1.f}, {0.5f,0.5f}};
    float[] vals = new float[]
    {  0.f,     1.f,      1.f,    0.5f };

    /****************
    double[][] pts = new double[][]
    {{0.,0.}, {1.,0}, {0.,1.}, {0.5,0.5}};
    double[] vals = new double[]
    {  0.,     1.,      1.,    0.5 };
    ****************/

    ScatterInterpSparse2 interp = new RbfScatterInterp(pts, vals, 0.2f);

    // walk along the diagonal.
    for( int i=0; i < 21; i++ ) {
      float frac = (float)i / 20.f;
      float x = frac;
      float y = frac;
      float v = interp.interp(x,y);
      System.out.println(v);
    }
    
  } //main

} //RbfScatterInterp
