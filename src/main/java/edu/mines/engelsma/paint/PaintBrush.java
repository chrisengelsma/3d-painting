/****************************************************************************
Copyright 2010, Colorado School of Mines and others.
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
****************************************************************************/
package edu.mines.engelsma.paint;

import edu.mines.engelsma.util.TimeSolver3;

import edu.mines.engelsma.vis.Contour;
import edu.mines.engelsma.vis.MarchingCubes;
import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.dsp.Tensors3;
import static edu.mines.jtk.util.MathPlus.max;
import static edu.mines.jtk.util.MathPlus.min;

/**
 * A 3D tensor-guided edu.mines.engelsma.paint brush.
 * This paintbrush paints voxels (3D pixels). Much like a 2D paintbrush, the
 * user must input the location and size of the brush. The 
 * edu.mines.engelsma.paint is then
 * guided by tensors by solving times using an eikonal equation.
 * @author Chris Engelsma and Dave Hale, Colorado School of Mines
 * @version 2010.12.10
 */
public class PaintBrush implements Painting3.DistanceMap {
  /**
   * Constructs a paintbrush.
   * @param n1 number of samples in 1st dimension (Z axis).
   * @param n2 number of samples in 2nd dimension (Y axis).
   * @param n3 number of samples in 3rd dimension (X axis).
   * @param pt painting tensors.
   */
  public PaintBrush(int n1, int n2, int n3, Tensors3 pt) {
    this(new Sampling(n1),new Sampling(n2),new Sampling(n3),pt);
  }

  /**
   * Constructs a paintbrush.
   * @param s1 sampling of 1st dimension (Z axis).
   * @param s2 sampling of 2nd dimension (Y axis).
   * @param s3 sampling of 3rd dimension (X axis).
   * @param pt painting tensors.
   */
  public PaintBrush(Sampling s1, Sampling s2, Sampling s3, Tensors3 pt) {
    _s1 = s1;
    _s2 = s2;
    _s3 = s3;
    _n1 = s1.getCount();
    _n2 = s2.getCount();
    _n3 = s3.getCount();
    _pt = pt;
    setSize(10);
  }

  /**
   * Sets the tensors for this paintbrush.
   * @param t painting tensors.
   */
  public void setTensors(Tensors3 t) {
    _pt = t;
  }

  /**
   * Gets the location of this brush.
   * @return location (k1,k2,k3) of this brush.
   */
  public int[] getLocation() {
    int[] k = new int[3];
    k[0] = _k1;
    k[1] = _k2;
    k[2] = _k3;
    return k;
  }

  /**
   * Sets the location of this brush.
   * The default location is (0,0,0).
   * @param k1 sample index in 1st dimension.
   * @param k2 sample index in 2nd dimension.
   * @param k3 sample index in 3rd dimension.
   */
  public void setLocation(int k1, int k2, int k3) {
    if (_k1!=k1 || _k2!=k2 || _k3!=k3) {
      _k1 = k1;
      _k2 = k2;
      _k3 = k3;
      _dirty = true;
    }
  }

  /**
   * Gets the radius size of the paintbrush.
   * @return the radius of the bounding sphere of this paintbrush.
   */
  public int getSize() {
    return _nh;
  }

  /**
   * Sets the size of this paintbrush.
   * For an identity painting tensor (for which time equals distance), the
   * brush size equals its radius. A brush with size zero covers only the one
   * sample where it is located. The default size is ten samples.
   * @param size the size.
   */
  public void setSize(int size) {
    _size = max(0,size);

    // Maximum time.
    _tmax = (float)(1+size);

    // Half of number of samples in array of brush times; modify as necessary.
    int nh = 1;
    while (nh<=size)
      nh += nh;

    // If number of samples in array of brush times has increased,
    // construct new brush tensors, time solver and marching cubes.
    if (_nh!=nh) {
      int nb = 1+2*nh;
      double db = 1.0;
      double fb = -nh;
      Sampling sb = new Sampling(nb,db,fb);
      _bt = new BrushTensors3();
      _ts = new TimeSolver3(nb,nb,nb,_bt);
      _mc = new MarchingCubes(sb,sb,sb,_ts.getTimes());
      _mc.setSwap13(true);
      _nb = nb;
      _nh = nh;
    }

    // Set maximum time for time solver, and note solution not valid.
    _ts.setMaxTime(2.0f*_tmax);
    _dirty = true;
  }

  /** 
   * Gets the contour for this paintbrush.
   * @return the contour.
   */
  public Contour getContour() {
    invokeTimeSolver();
    double d1 = _s1.getDelta(), d2 = _s2.getDelta(), d3 = _s3.getDelta();
    double f1 = _s1.getFirst(), f2 = _s2.getFirst(), f3 = _s3.getFirst();
    
    //compute on subset
    Sampling s1 = new Sampling(_nb,d1,f1+(_k1-_nh)*d1);
    Sampling s2 = new Sampling(_nb,d2,f2+(_k2-_nh)*d2);
    Sampling s3 = new Sampling(_nb,d3,f3+(_k3-_nh)*d3);
    _mc.setSampling(s1,s2,s3);
    
    Contour contour = _mc.getContour((float)_size);
    _contour = new Contour();
    _contour.x = contour.x; // contour vertices
    _contour.u = contour.u; // contour normals
    _contour.i = contour.i; // contour vertex indices
    return _contour;
  }

  /**
   * Gets the distance from any point (i1,i2,i3) to the origin of the brush.
   * For an identity painting tensor, this is Euclidean distance.
   * @param i1 sample index in 1st dimension.
   * @param i2 sample index in 2nd dimension.
   * @param i3 sample index in 3rd dimension.
   * @return the distance from (i1,i2,i3) to the brush's origin.
   */
  public float getDistance(int i1, int i2, int i3) {
    int ii1 = i1-_k1+_nh;
    int ii2 = i2-_k2+_nh;
    int ii3 = i3-_k3+_nh;

    if (ii1<0 || ii1>=_nb) return INFINITY;
    if (ii2<0 || ii2>=_nb) return INFINITY;
    if (ii3<0 || ii3>=_nb) return INFINITY;

    invokeTimeSolver();

    float time = _ts.getTimes()[ii3][ii2][ii1];
    return time;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  // Brush tensors are a subset of the painting tensors.
  private class BrushTensors3 implements Tensors3 {
    public void getTensor(int i1, int i2, int i3, float[] a) {
      i1 = max(0,min(_n1-1,i1+_k1-_nh));
      i2 = max(0,min(_n2-1,i2+_k2-_nh));
      i3 = max(0,min(_n3-1,i3+_k3-_nh));
      _pt.getTensor(i1,i2,i3,a);
    }
  }

  private int _n1,_n2,_n3;
  private Sampling _s1,_s2,_s3;
  private Tensors3 _pt;
  private int _k1,_k2,_k3;
  private int _size,_nb,_nh;
  private float _tmax;
  private MarchingCubes _mc;
  private BrushTensors3 _bt;
  private TimeSolver3 _ts;
  private boolean _dirty;
  private Contour _contour;
  private static final float INFINITY = Float.MAX_VALUE;

  // Debugging tool
  private void trace(Object o) {
    System.out.println(o);
  }

  // Called when times need to be recomputed
  private void invokeTimeSolver() {
    if (_dirty) {
      _ts.reset();
      _ts.zeroAt(_nh,_nh,_nh);
      _dirty = false;
    }
  }
}
