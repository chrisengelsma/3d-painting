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

import edu.mines.engelsma.vis.Contour;
import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.util.Parallel;
import static edu.mines.jtk.util.ArrayMath.*;

import java.io.Serializable;
import java.util.HashMap;

/**
 * A 3D painting.
 * 3D painting is intended to be an intuitive extension of 2D paintings in 3D.
 * A 3D painting can be thought of as a volume of painted voxels (3D pixels).
 * <p>
 * The painting data structure maintains intersections between grid points,
 * representing each 3D painting with sub-voxel precision. These intersection
 * points are stored as bytes [0,100) in three separate arrays. For example,
 * v1[0][0][0] stores the edge intersection between (0,0,0) and (0,0,1),
 * v2[0][0][0] stores the edge intersection between (0,0,0) and (0,1,0),
 * and v3[0][0][0] stores the edge intersection between (0,0,0) and (1,0,0).
 * If no intersection exists, -1 is assigned.
 * @author Chris Engelsma, Colorado School of Mines
 * @version 2010.12.10
 */
public class Painting3 implements Serializable {
  private static final long serialVersionUID = 1L;

  /**
   * The concurrency of this process.
   */
  public enum Concurrency {
    PARALLEL,
    SERIAL
  };
  
  /**
   * Constructs a new 3D Painting.
   * @param n3 size of image in dimension 3.
   * @param n2 size of image in dimension 2.
   * @param n1 size of image in dimension 1.
   */
  public Painting3(int n1, int n2, int n3) {
    this(new float[n3][n2][n1]);
  }

  /**
   * Constructs a new 3D Painting.
   * @param paint a 3D image.
   */
  public Painting3(float[][][] paint) {
    this(new Sampling(paint[0][0].length),
         new Sampling(paint[0].length),
         new Sampling(paint.length),
         paint);
  }

  /**
   * Constructs a new 3D Painting.
   * @param s1 sampling in the 1st dimension.
   * @param s2 sampling in the 2nd dimension.
   * @param s3 sampling in the 3rd dimension.
   */
  public Painting3(
    Sampling s1, Sampling s2, Sampling s3)
  {
    this(s1,s2,s3,new float[s3.getCount()][s2.getCount()][s1.getCount()]);
  }
  
  /**
   * Constructs a new 3D Painting.
   * @param s1 sampling in the 1st dimension.
   * @param s2 sampling in the 2nd dimension.
   * @param s3 sampling in the 3rd dimension.
   * @param paint a 3D image.
   */
  public Painting3(
    Sampling s1, Sampling s2, Sampling s3, float[][][] paint)
  {
    _paint = paint;

    _n3 = s3.getCount();
    _n2 = s2.getCount();
    _n1 = s1.getCount();

    _f3 = s3.getFirst();
    _f2 = s2.getFirst();
    _f1 = s1.getFirst();

    _d3 = s3.getDelta();
    _d2 = s2.getDelta();
    _d1 = s1.getDelta();

    _v3 = new byte[_n3][_n2][_n1];
    _v2 = new byte[_n3][_n2][_n1];
    _v1 = new byte[_n3][_n2][_n1];

    initializeAllArrays(false);
  }

  /**
   * Paints the image at a given coordinate.
   * @param i1 the index in the 1st dimension.
   * @param i2 the index in the 2nd dimension.
   * @param i3 the index in the 3rd dimension.
   * @param v the value to edu.mines.engelsma.paint.
   * @param d the maximum size.
   * @param dm the distance map.
   */
  public void paintAt(
    int i1, int i2, int i3,
    float v, float d, Painting3.DistanceMap dm)
  {
    dm.setLocation(i1,i2,i3);
    markBorders(i3,i2,i1,v,d,dm,false);
  }

  /**
   * Erases edu.mines.engelsma.paint from the image at a given coordinate.
   * @param i1 the index in the 1st dimension. 
   * @param i2 the index in the 2nd dimension.
   * @param i3 the index in the 3rd dimension.
   * @param d the maximum size.
   * @param dm the distance map.
   */
  public void eraseAt(
    int i1, int i2, int i3,
    float d, Painting3.DistanceMap dm)
  {
    dm.setLocation(i1,i2,i3);
    markBorders(i3,i2,i1,0,d,dm,true);
  }

  /**
   * Sets all painted voxels to 0.
   */
  public void eraseAll() {
    initializeAllArrays(true);
  }

  /**
   * Returns the painted image.
   * @return the painted image.
   */
  public float[][][] getPaint() {
    return _paint;
  }

  /**
   * Swaps the first and third index.
   * @param swap true, if swapping indices 1 and 3; false, otherwise.
   */
  public void setSwap13(boolean swap) {
    _swap13 = swap;
  }

  /**
   * Sets the concurrency of this procedure.
   * @param c the concurrency.
   */
  public void setConcurrency(Concurrency c) {
    _concurrency = c;
  }

  /**
   * Gets the 3D contour of painted voxels.
   * @param c the contour value.
   * @return the 3D contour of painted voxels.
   */
  public Contour getContour(float c) {

    // Run marching cubes.
    initiateMarch(c);

    int xlen = 0,tlen = 0,stlen = 0;
    for (int i=0; i<_n3; ++i) {
      if (_xlist[i]!=null) xlen+=_xlist[i].n;
      if (_tlist[i]!=null) tlen+=_tlist[i].n;
      if (_stlist[i]!=null) stlen+=_stlist[i].n;
    }

    /* Because the marching cubes algorithm can be run in parallel, 
     * the resulting list of vertices needs to be concatenated. This is
     * currently run in serial, and is a bottleneck. Ideally, this should
     * be parallelized.
     */
    Contour contour = new Contour();
    contour.x = new float[xlen];
    contour.i = new int[tlen];
    int[] slist = new int[stlen];
    int xn = 0,tn = 0,stn = 0;
    for (int i=0; i<_n3; ++i) {
      if (_xlist[i]!=null) {
        System.arraycopy(_xlist[i].trim(),0,contour.x,xn,_xlist[i].n);
        xn+=_xlist[i].n;
      }
      if (_tlist[i]!=null) {
        System.arraycopy(_tlist[i].trim(),0,contour.i,tn,_tlist[i].n);
        tn+=_tlist[i].n;
      }
      if (_stlist[i]!=null) {
        int sn = _stlist[i].n;
        System.arraycopy(_stlist[i].trim(),0,slist,stn,sn);
        stn+=sn;
      }
    }

    HashMap<Integer,Integer> hm = new HashMap<Integer,Integer>();
    int key,value = 0;
    for (int i=0; i<stn; ++i) {
      key = slist[i];
      if (!hm.containsKey(key)) {
        hm.put(key,value++);
      }
    }
    
    for (int i=0; i<tlen; ++i) {
      key = contour.i[i];
      contour.i[i] = hm.get(key);
    }

    // Swap the 1st and 3rd dimension if necessary.
    if (_swap13) {
      float[] x = contour.x;
      float[] u = contour.u;
      for (int i=x.length-3; i>=0; i-=3) {
        float x1 = x[i  ];
        float x3 = x[i+2];
        x[i  ] = x3;
        x[i+2] = x1;
      }
    }
    return contour;
  }

  /**
   * Gets the edge intersections at a given sample.
   * Edge intersections are returned as bytes in the range [-1,100), where
   * a -1 value means that no intersection exist.
   * @param i1 sample in the 1st dimension.
   * @param i2 sample in the 2nd dimension.
   * @param i3 sample in the 3rd dimension.
   * @return the edge intersections (e1,e2,e3).
   */
  public byte[] getEdgeIntersectionsAt(int i1, int i2, int i3) {
    if (i1<0) i1 = 0; if (i1>_n1-1) i1 = _n1-1;
    if (i2<0) i2 = 0; if (i2>_n2-1) i2 = _n2-1;
    if (i3<0) i3 = 0; if (i3>_n3-1) i3 = _n3-1;
    byte[] e = new byte[3];
    e[0] = _v1[i3][i2][i1];
    e[1] = _v2[i3][i2][i1];
    e[2] = _v3[i3][i2][i1];
    return e;
  }

  /**
   * Sets the edge intersection in the first dimension at a given sample.
   * Edge intersections will be clipped between the range [-1,100).
   * @param i1 sample in the 1st dimension.
   * @param i2 sample in the 2nd dimension.
   * @param i3 sample in the 3rd dimension.
   * @param v a byte value within the range [-1,100).
   */
  public void setEdgeIntersection1At(int i1, int i2, int i3, byte v) {
    v = (byte)max(-1,min(100,v));
    _v1[i3][i2][i1] = v;
  }

  /**
   * Sets the edge intersection in the second dimension at a given sample.
   * Edge intersections will be clipped between the range [-1,100).
   * @param i1 sample in the 1st dimension.
   * @param i2 sample in the 2nd dimension.
   * @param i3 sample in the 3rd dimension.
   * @param v a byte value within the range [-1,100).
   */
  public void setEdgeIntersection2At(int i1, int i2, int i3, byte v) {
    v = (byte)max(-1,min(100,v));
    _v2[i3][i2][i1] = v;
  }

  /**
   * Sets the edge intersection in the third dimension at a given sample.
   * Edge intersections will be clipped between the range [-1,100).
   * @param i1 sample in the 1st dimension.
   * @param i2 sample in the 2nd dimension.
   * @param i3 sample in the 3rd dimension.
   * @param v a byte value within the range [-1,100).
   */
  public void setEdgeIntersection3At(int i1, int i2, int i3, byte v) {
    v = (byte)max(-1,min(100,v));
    _v3[i3][i2][i1] = v;
  }
  
  ///////////////////////////////////////////////////////////////////////////
  // private

  private class IntList {
    public int n;      // the number of submitted values.
    public int[] a = new int[32]; // the array of ints.
    public void add(int i) {
      if (n==a.length) {
        int[] t = new int[2*n];
        System.arraycopy(a,0,t,0,n);
        a = t;
      }
      a[n++] = i;
    }
    public int[] trim() {
      int[] t = new int[n];
      System.arraycopy(a,0,t,0,n);
      return t;
    }
  }

  private class FloatList {
    public int n;      // the number of submitted values.
    public float[] a = new float[32]; // the array of floats.
    public void add(float d) {
      if (n==a.length) {
        float[] t = new float[2*n];
        System.arraycopy(a,0,t,0,n);
        a = t;
      }
      a[n++] = d;
    }
    public float[] trim() {
      float[] t = new float[n];
      System.arraycopy(a,0,t,0,n);
      return t;
    }
  }
  
  private int _n3,_n2,_n1;
  private double _f3,_f2,_f1;
  private double _d3,_d2,_d1;

  private FloatList[] _xlist; // List of vertices.
  private IntList[] _tlist;   // List of triangle indices.
  private IntList[] _stlist;  // List of ordered triangle indices.
  
  private float[][][] _paint;
  private int[][][] _ix0;
  private int[][][] _ix1;
  private int[][][] _ix2;
  private byte[][][] _v3; // edge in dimension 3
  private byte[][][] _v2; // edge in dimension 2
  private byte[][][] _v1; // edge in dimension 1

  private boolean _swap13 = true; // Swap the 1st and 3rd dimensions?

//  private Concurrency _concurrency = Concurrency.PARALLEL;
  private Concurrency _concurrency = Concurrency.PARALLEL;

  public interface DistanceMap {
    float getDistance(int i1, int i2, int i3);
    void setLocation(int i1, int i2, int i3);
    int getSize();
  }

  /**
   * Marks the borders of the painted area.
   */
  private void markBorders(
    int ii3, int ii2, int ii1,
    float v, float dmax, Painting3.DistanceMap dm,
    boolean erasing)
  {
    if (ii3<0 || ii2<0 || ii1<0) return;
    if (ii3>_n3-1 || ii2>_n2-1 || ii1>_n1-1) return;

    int nh = dm.getSize();

    int k1 = max(0,min(_n1,ii1-nh));
    int k2 = max(0,min(_n2,ii2-nh));
    int k3 = max(0,min(_n3,ii3-nh));

    int j1 = max(0,min(_n1,ii1+nh));
    int j2 = max(0,min(_n2,ii2+nh));
    int j3 = max(0,min(_n3,ii3+nh));

    // Compute on a subcube
    for (int i3=k3; i3<j3; ++i3) {
      for (int i2=k2; i2<j2; ++i2) {
        for (int i1=k1; i1<j1; ++i1) {

          /* Get distances for this sample and the succeeding samples in the
           * 1st, 2nd and 3rd dimensions.
           */
          float d0 = dm.getDistance(i1,i2,i3);
          float d1 = (i1!=_n1-1)?dm.getDistance(i1+1,i2  ,i3  ):d0;
          float d2 = (i2!=_n2-1)?dm.getDistance(i1  ,i2+1,i3  ):d0;
          float d3 = (i3!=_n3-1)?dm.getDistance(i1  ,i2  ,i3+1):d0;

          /* If this sample is within the painted area, 
          edu.mines.engelsma.paint it. */
          if (d0<=dmax) {
            if(!erasing) _paint[i3][i2][i1] = v;
            else         _paint[i3][i2][i1] = 0.0f;
          }

          /* Find edges that intersect and compute that intersection's
           * percent-distance [0,100] using linear interpolation.
           *
           * Conditions are as follows for edge assignment:
           * Assuming we have an edge, we have four possible scenarios:
           *
           * o----------------o : No painted value.
           *                      Do nothing.
           *
           * o-----x----------o : First time being painted.
           *                      Assign the edge value.
           *
           * o-----x----X-----o : Repaint, shifting value up.
           *       ----->         Check if the seed is painted.
           *                      If true, accept the new value.
           *
           * o-----X----x-----o : Repaint, shifting value down.
           *       <-----         Check if the seed is painted.
           *                      If true, do nothing.
           */
          byte v0,nv;
          float p  = _paint[i3  ][i2  ][i1  ];

          if (i1!=_n1-1) {
            float p1 = _paint[i3  ][i2  ][i1+1];
            if (intersectsEdge(dmax,d0,d1)) {
              v0 = _v1[i3][i2][i1];
              nv = (byte)min(100.0f,max(0.0f,(dmax-d0)/(d1-d0)*100));
              if (v0==-1) _v1[i3][i2][i1] = nv;
              else {
                if ((nv>v0 && p!=0.0f) || (nv<v0 && p==0.0f))
                  _v1[i3][i2][i1] = nv;
              }
            } else {
              if (p==p1) _v1[i3][i2][i1] = -1;
            }
          }

          if (i2!=_n2-1) {
            float p2 = _paint[i3  ][i2+1][i1  ];
            if (intersectsEdge(dmax,d0,d2)) {
              v0 = _v2[i3][i2][i1];
              nv = (byte)min(100.0f,max(0.0f,(dmax-d0)/(d2-d0)*100));
              if (v0==-1) _v2[i3][i2][i1] = nv;
              else {
                if ((nv>v0 && p!=0.0f) || (nv<v0 && p==0.0f))
                  _v2[i3][i2][i1] = nv;
              }
            } else {
              if (p==p2) _v2[i3][i2][i1] = -1;
            }
          }

          if (i3!=_n3-1) {
            float p3 = _paint[i3+1][i2  ][i1  ];
            if (intersectsEdge(dmax,d0,d3)) {
              v0 = _v3[i3][i2][i1];
              nv = (byte)min(100.0f,max(0.0f,(dmax-d0)/(d3-d0)*100));
              if (v0==-1) _v3[i3][i2][i1] = nv;
              else {
                if ((nv>v0 && p!=0.0f) || (nv<v0 && p==0.0f))
                  _v3[i3][i2][i1] = nv;
              }
            } else {
              if (p==p3) _v3[i3][i2][i1] = -1;
            }
          }
        }
      }
    }
  }

  /**
   * Determines whether an edge is intersected from this sample.
   */
  private boolean intersectsEdge(float dmax, float d0, float d1) {
    if (d0-dmax<0 && d1-dmax>=0) return true;
    if (d0-dmax>0 && d1-dmax<=0) return true;
    return false;
  }

  /**
   * Extracts a painted contour value in order to render in 3D.
   */
  private float[][][] isolateContourValue(float c) {
    final float val = c;
    final float[][][] array = new float[_n3][_n2][_n1];
    Parallel.loop(_n3-1,new Parallel.LoopInt() {
      public void compute(int i3) {
        isolateSlab(i3,val,array);
      }
    });

    return array;
  }

  /**
   * Extracts a single slab (used for parallel concurrency).
   */
  private void isolateSlab(int i3, float val, float[][][] array) {
    for (int i2=0; i2<_n2-1; ++i2)
      for (int i1=0; i1<_n1-1; ++i1)
        array[i3][i2][i1] = (_paint[i3][i2][i1]==val)?val:0.0f;
  }

  /**
   * Sets up for the marching cubes algorithm.
   */
  private void initiateMarch(float c) {
    float[][][] array = isolateContourValue(c);
    _xlist = new FloatList[_n3];
    _tlist = new IntList[_n3];
    _stlist = new IntList[_n3];

    _ix0 = new int[_n3][_n2][_n1];
    _ix1 = new int[_n3][_n2][_n1];
    _ix2 = new int[_n3][_n2][_n1];
    for (int i3=0; i3<_n3; ++i3) {
      for (int i2=0; i2<_n2; ++i2) {
        for (int i1=0; i1<_n1; ++i1) {
          _ix0[i3][i2][i1] = -1;
          _ix1[i3][i2][i1] = -1;
          _ix2[i3][i2][i1] = -1;
        }
      }
    }

    if (_concurrency==Concurrency.PARALLEL) marchParallel(array,c);
    else marchSerial(array,c);

  }

  /**
   * Performs an adapted version of marching cubes whereby the precomputed
   * edge-intersections are employed.
   * Serial version.
   */
  private void marchSerial(float[][][] array, float c) {
    int i3step = 2;
    int i3start = 0;
    int i3stop = _n3-1;
    for (int i3pass=0; i3pass<i3step; ++i3pass,++i3start) {
      for (int i3=i3start; i3<i3stop; i3+=i3step) {
        _xlist[i3] = new FloatList();
        _tlist[i3] = new IntList();
        _stlist[i3] = new IntList();
        march(array,c,i3,_xlist[i3],_tlist[i3],_stlist[i3]);
      }
    }
  }

  /**
   * Performs an adapted version of marching cubes whereby the precomputed
   * edge-intersections are employed.
   * Parallel version.
   */
  private void marchParallel(float[][][] array, float c) {
    final float[][][] farray = array;
    final float fc = c;
    final int i3step = 2;
    int i3start = 0;
    final int i3stop = _n3-1;
    for (int i3pass=0; i3pass<i3step; ++i3pass,++i3start) {
      Parallel.loop(i3start,i3stop,i3step,new Parallel.LoopInt() {
        public void compute(int i3) {
          _tlist[i3] = new IntList();
          _stlist[i3] = new IntList();
          _xlist[i3] = new FloatList();
          march(farray,fc,i3,_xlist[i3],_tlist[i3],_stlist[i3]);
        }
      });
    }
  }

  /**
   * Produces a packed array of vertices and indices using marching cubes
   * logic. Also, precomputed edge intersections are used.
   */
  private void march(
    float[][][] paint, float c, int i3,
    FloatList xlist, IntList tlist, IntList stlist)
  {
    int nx3 = _n1*_n2*i3;

    for (int i2=0; i2<_n2-1; ++i2) {
      for (int i1=0; i1<_n1-1; ++i1) {

        // Eight corner values for this cube.
        float c0 = paint[i3  ][i2  ][i1  ];
        float c1 = paint[i3  ][i2  ][i1+1];
        float c2 = paint[i3  ][i2+1][i1+1];
        float c3 = paint[i3  ][i2+1][i1  ];
        float c4 = paint[i3+1][i2  ][i1  ];
        float c5 = paint[i3+1][i2  ][i1+1];
        float c6 = paint[i3+1][i2+1][i1+1];
        float c7 = paint[i3+1][i2+1][i1  ];

        // Case index
        int ci = 0;
        if (c0==c) ci +=   1;
        if (c1==c) ci +=   2;
        if (c2==c) ci +=   4;
        if (c3==c) ci +=   8;
        if (c4==c) ci +=  16;
        if (c5==c) ci +=  32;
        if (c6==c) ci +=  64;
        if (c7==c) ci += 128;

        // If at least one triangle for this case, ...
        if (ci>0 && ci<255) {

          // Edges intersected by contour.
          int[] edges = _edges[ci];
          int ne = edges.length;

          // For each edge triplet
          for (int ie=0; ie<ne; ie+=3) {

            // For each edge
            for (int je=0; je<3; ++je) {
              // Pull the interpolated value from the byte arrays
              int edge = edges[ie+je];
              float p;
              int j1,j2,j3,kk;
              switch(edge) {
              case 0: // 0->1
                j1 = i1;
                j2 = i2;
                j3 = i3;
                kk = 0;
                break;
              case 1: // 1->2
                j1 = i1+1;
                j2 = i2;
                j3 = i3;
                kk = 1;
                break;
              case 2: // 3->2
                j1 = i1;
                j2 = i2+1;
                j3 = i3;
                kk = 0;
                break;
              case 3: // 0->3
                j1 = i1;
                j2 = i2;
                j3 = i3;
                kk = 1;
                break;
              case 4: // 4->5
                j1 = i1;
                j2 = i2;
                j3 = i3+1;
                kk = 0;
                break;
              case 5: // 5->6
                j1 = i1+1;
                j2 = i2;
                j3 = i3+1;
                kk = 1;
                break;
              case 6: // 7->6
                j1 = i1;
                j2 = i2+1;
                j3 = i3+1;
                kk = 0;
                break;
              case 7: // 4->7
                j1 = i1;
                j2 = i2;
                j3 = i3+1;
                kk = 1;
                break;
              case 8: // 0->4
                j1 = i1;
                j2 = i2;
                j3 = i3;
                kk = 2;
                break;
              case 9: // 1->5
                j1 = i1+1;
                j2 = i2;
                j3 = i3;
                kk = 2;
                break;
              case 10: // 3->7
                j1 = i1;
                j2 = i2+1;
                j3 = i3;
                kk = 2;
                break;
              default: // 2->6
                j1 = i1+1;
                j2 = i2+1;
                j3 = i3;
                kk = 2;
              }

              /*
               * Pre-computed edge intersections are used here.
               * Depending on the dimension (kk), the proper edge intersection
               * is pulled from its corresponding array.
               */
              if (kk==0)      p = (float)_v1[j3][j2][j1]/100.0f;
              else if (kk==1) p = (float)_v2[j3][j2][j1]/100.0f;
              else            p = (float)_v3[j3][j2][j1]/100.0f;

              // Index of vertex, if already computed; or -1, if not yet.
              int ix;
              if (kk==0)      ix = _ix0[j3][j2][j1];
              else if (kk==1) ix = _ix1[j3][j2][j1];
              else            ix = _ix2[j3][j2][j1];
              
              // If vertex not yet stored...
              if (ix<0) {
                double x1,x2,x3;
                switch(kk) {
                case 0:  // axis 1
                  x1 = _f1+_d1*(j1+p);
                  x2 = _f2+_d2*(j2  );
                  x3 = _f3+_d3*(j3  );
                  break;
                case 1:  // axis 2
                  x1 = _f1+_d1*(j1  );
                  x2 = _f2+_d2*(j2+p);
                  x3 = _f3+_d3*(j3  );
                  break;
                default: // axis 3
                  x1 = _f1+_d1*(j1  );
                  x2 = _f2+_d2*(j2  );
                  x3 = _f3+_d3*(j3+p);
                }
                ix = nx3++;
                xlist.add((float)x1);
                xlist.add((float)x2);
                xlist.add((float)x3);
                stlist.add(ix);
                if (kk==0)      _ix0[j3][j2][j1] = ix;
                else if (kk==1) _ix1[j3][j2][j1] = ix;
                else            _ix2[j3][j2][j1] = ix;
              }
              tlist.add(ix);
            }
          }
        }
      }
    }
  }

  /**
   * Resets all arrays. Paint array is set to zeros, edge intersection arrays
   * are set to negative ones.
   */
  private void initializeAllArrays(boolean all) {
    for (int i3=0; i3<_n3; ++i3) {
      for (int i2=0; i2<_n2; ++i2) {
        for (int i1=0; i1<_n1; ++i1) {
          if (all) _paint[i3][i2][i1] = 0.0f;
          _v1[i3][i2][i1] = -1;
          _v2[i3][i2][i1] = -1;
          _v3[i3][i2][i1] = -1;
        }
      }
    }
  }

  /**
   * Debugging.
   */
  private static void print(Object o) {
    System.out.println(o);
  }

  /* Edges intersected. Each group of three indices represents a triangle.
   * For the eight sample values in each cube, there are 256 cases. However,
   * most of those 256 cases are complements or rotations of 16 base cases.
   * Comments at end of each line are case number and base-case number.
   * This table was adopted from one in VTK, the Visualization Toolkit.
   */
  private static final int[][] _edges = {
    {}, // 0 0
    { 0, 3, 8}, // 1 1
    { 0, 9, 1}, // 2 1
    { 1, 3, 8, 9, 1, 8}, // 3 2
    { 1, 11, 2}, // 4 1
    { 0, 3, 8, 1, 11, 2}, // 5 3
    { 9, 11, 2, 0, 9, 2}, // 6 2
    { 2, 3, 8, 2, 8, 11, 11, 8, 9}, // 7 5
    { 3, 2, 10}, // 8 1
    { 0, 2, 10, 8, 0, 10}, // 9 2
    { 1, 0, 9, 2, 10, 3}, // 10 3
    { 1, 2, 10, 1, 10, 9, 9, 10, 8}, // 11 5
    { 3, 1, 11, 10, 3, 11}, // 12 2
    { 0, 1, 11, 0, 11, 8, 8, 11, 10}, // 13 5
    { 3, 0, 9, 3, 9, 10, 10, 9, 11}, // 14 5
    { 9, 11, 8, 11, 10, 8}, // 15 8
    { 4, 8, 7}, // 16 1
    { 4, 0, 3, 7, 4, 3}, // 17 2
    { 0, 9, 1, 8, 7, 4}, // 18 3
    { 4, 9, 1, 4, 1, 7, 7, 1, 3}, // 19 5
    { 1, 11, 2, 8, 7, 4}, // 20 4
    { 3, 7, 4, 3, 4, 0, 1, 11, 2}, // 21 7
    { 9, 11, 2, 9, 2, 0, 8, 7, 4}, // 22 7
    { 2, 9, 11, 2, 7, 9, 2, 3, 7, 7, 4, 9}, // 23 14
    { 8, 7, 4, 3, 2, 10}, // 24 3
    {10, 7, 4, 10, 4, 2, 2, 4, 0}, // 25 5
    { 9, 1, 0, 8, 7, 4, 2, 10, 3}, // 26 6
    { 4, 10, 7, 9, 10, 4, 9, 2, 10, 9, 1, 2}, // 27 9
    { 3, 1, 11, 3, 11, 10, 7, 4, 8}, // 28 7
    { 1, 11, 10, 1, 10, 4, 1, 4, 0, 7, 4, 10}, // 29 11
    { 4, 8, 7, 9, 10, 0, 9, 11, 10, 10, 3, 0}, // 30 12
    { 4, 10, 7, 4, 9, 10, 9, 11, 10}, // 31 5
    { 9, 4, 5}, // 32 1
    { 9, 4, 5, 0, 3, 8}, // 33 3
    { 0, 4, 5, 1, 0, 5}, // 34 2
    { 8, 4, 5, 8, 5, 3, 3, 5, 1}, // 35 5
    { 1, 11, 2, 9, 4, 5}, // 36 3
    { 3, 8, 0, 1, 11, 2, 4, 5, 9}, // 37 6
    { 5, 11, 2, 5, 2, 4, 4, 2, 0}, // 38 5
    { 2, 5, 11, 3, 5, 2, 3, 4, 5, 3, 8, 4}, // 39 9
    { 9, 4, 5, 2, 10, 3}, // 40 4
    { 0, 2, 10, 0, 10, 8, 4, 5, 9}, // 41 7
    { 0, 4, 5, 0, 5, 1, 2, 10, 3}, // 42 7
    { 2, 5, 1, 2, 8, 5, 2, 10, 8, 4, 5, 8}, // 43 11
    {11, 10, 3, 11, 3, 1, 9, 4, 5}, // 44 7
    { 4, 5, 9, 0, 1, 8, 8, 1, 11, 8, 11, 10}, // 45 12
    { 5, 0, 4, 5, 10, 0, 5, 11, 10, 10, 3, 0}, // 46 14
    { 5, 8, 4, 5, 11, 8, 11, 10, 8}, // 47 5
    { 9, 8, 7, 5, 9, 7}, // 48 2
    { 9, 0, 3, 9, 3, 5, 5, 3, 7}, // 49 5
    { 0, 8, 7, 0, 7, 1, 1, 7, 5}, // 50 5
    { 1, 3, 5, 3, 7, 5}, // 51 8
    { 9, 8, 7, 9, 7, 5, 11, 2, 1}, // 52 7
    {11, 2, 1, 9, 0, 5, 5, 0, 3, 5, 3, 7}, // 53 12
    { 8, 2, 0, 8, 5, 2, 8, 7, 5, 11, 2, 5}, // 54 11
    { 2, 5, 11, 2, 3, 5, 3, 7, 5}, // 55 5
    { 7, 5, 9, 7, 9, 8, 3, 2, 10}, // 56 7
    { 9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 10, 7}, // 57 14
    { 2, 10, 3, 0, 8, 1, 1, 8, 7, 1, 7, 5}, // 58 12
    {10, 1, 2, 10, 7, 1, 7, 5, 1}, // 59 5
    { 9, 8, 5, 8, 7, 5, 11, 3, 1, 11, 10, 3}, // 60 10
    { 5, 0, 7, 5, 9, 0, 7, 0, 10, 1, 11, 0, 10, 0, 11}, // 61 7
    {10, 0, 11, 10, 3, 0, 11, 0, 5, 8, 7, 0, 5, 0, 7}, // 62 7
    {10, 5, 11, 7, 5, 10}, // 63 2
    {11, 5, 6}, // 64 1
    { 0, 3, 8, 5, 6, 11}, // 65 4
    { 9, 1, 0, 5, 6, 11}, // 66 3
    { 1, 3, 8, 1, 8, 9, 5, 6, 11}, // 67 7
    { 1, 5, 6, 2, 1, 6}, // 68 2
    { 1, 5, 6, 1, 6, 2, 3, 8, 0}, // 69 7
    { 9, 5, 6, 9, 6, 0, 0, 6, 2}, // 70 5
    { 5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2}, // 71 11
    { 2, 10, 3, 11, 5, 6}, // 72 3
    {10, 8, 0, 10, 0, 2, 11, 5, 6}, // 73 7
    { 0, 9, 1, 2, 10, 3, 5, 6, 11}, // 74 6
    { 5, 6, 11, 1, 2, 9, 9, 2, 10, 9, 10, 8}, // 75 12
    { 6, 10, 3, 6, 3, 5, 5, 3, 1}, // 76 5
    { 0, 10, 8, 0, 5, 10, 0, 1, 5, 5, 6, 10}, // 77 14
    { 3, 6, 10, 0, 6, 3, 0, 5, 6, 0, 9, 5}, // 78 9
    { 6, 9, 5, 6, 10, 9, 10, 8, 9}, // 79 5
    { 5, 6, 11, 4, 8, 7}, // 80 3
    { 4, 0, 3, 4, 3, 7, 6, 11, 5}, // 81 7
    { 1, 0, 9, 5, 6, 11, 8, 7, 4}, // 82 6
    {11, 5, 6, 1, 7, 9, 1, 3, 7, 7, 4, 9}, // 83 12
    { 6, 2, 1, 6, 1, 5, 4, 8, 7}, // 84 7
    { 1, 5, 2, 5, 6, 2, 3, 4, 0, 3, 7, 4}, // 85 10
    { 8, 7, 4, 9, 5, 0, 0, 5, 6, 0, 6, 2}, // 86 12
    { 7, 9, 3, 7, 4, 9, 3, 9, 2, 5, 6, 9, 2, 9, 6}, // 87 7
    { 3, 2, 10, 7, 4, 8, 11, 5, 6}, // 88 6
    { 5, 6, 11, 4, 2, 7, 4, 0, 2, 2, 10, 7}, // 89 12
    { 0, 9, 1, 4, 8, 7, 2, 10, 3, 5, 6, 11}, // 90 13
    { 9, 1, 2, 9, 2, 10, 9, 10, 4, 7, 4, 10, 5, 6, 11}, // 91 6
    { 8, 7, 4, 3, 5, 10, 3, 1, 5, 5, 6, 10}, // 92 12
    { 5, 10, 1, 5, 6, 10, 1, 10, 0, 7, 4, 10, 0, 10, 4}, // 93 7
    { 0, 9, 5, 0, 5, 6, 0, 6, 3, 10, 3, 6, 8, 7, 4}, // 94 6
    { 6, 9, 5, 6, 10, 9, 4, 9, 7, 7, 9, 10}, // 95 3
    {11, 9, 4, 6, 11, 4}, // 96 2
    { 4, 6, 11, 4, 11, 9, 0, 3, 8}, ///97 7
    {11, 1, 0, 11, 0, 6, 6, 0, 4}, // 98 5
    { 8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 11, 1}, // 99 14
    { 1, 9, 4, 1, 4, 2, 2, 4, 6}, // 100 5
    { 3, 8, 0, 1, 9, 2, 2, 9, 4, 2, 4, 6}, // 101 12
    { 0, 4, 2, 4, 6, 2}, // 102 8
    { 8, 2, 3, 8, 4, 2, 4, 6, 2}, // 103 5
    {11, 9, 4, 11, 4, 6, 10, 3, 2}, // 104 7
    { 0, 2, 8, 2, 10, 8, 4, 11, 9, 4, 6, 11}, // 105 10
    { 3, 2, 10, 0, 6, 1, 0, 4, 6, 6, 11, 1}, // 106 12
    { 6, 1, 4, 6, 11, 1, 4, 1, 8, 2, 10, 1, 8, 1, 10}, // 107 7
    { 9, 4, 6, 9, 6, 3, 9, 3, 1, 10, 3, 6}, // 108 11
    { 8, 1, 10, 8, 0, 1, 10, 1, 6, 9, 4, 1, 6, 1, 4}, // 109 7
    { 3, 6, 10, 3, 0, 6, 0, 4, 6}, // 110 5
    { 6, 8, 4, 10, 8, 6}, // 111 2
    { 7, 6, 11, 7, 11, 8, 8, 11, 9}, // 112 5
    { 0, 3, 7, 0, 7, 11, 0, 11, 9, 6, 11, 7}, // 113 11
    {11, 7, 6, 1, 7, 11, 1, 8, 7, 1, 0, 8}, // 114 9
    {11, 7, 6, 11, 1, 7, 1, 3, 7}, // 115 5
    { 1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6}, // 116 14
    { 2, 9, 6, 2, 1, 9, 6, 9, 7, 0, 3, 9, 7, 9, 3}, // 117 7
    { 7, 0, 8, 7, 6, 0, 6, 2, 0}, // 118 5
    { 7, 2, 3, 6, 2, 7}, // 119 2
    { 2, 10, 3, 11, 8, 6, 11, 9, 8, 8, 7, 6}, // 120 12
    { 2, 7, 0, 2, 10, 7, 0, 7, 9, 6, 11, 7, 9, 7, 11}, // 121 7
    { 1, 0, 8, 1, 8, 7, 1, 7, 11, 6, 11, 7, 2, 10, 3}, // 122 6
    {10, 1, 2, 10, 7, 1, 11, 1, 6, 6, 1, 7}, // 123 3
    { 8, 6, 9, 8, 7, 6, 9, 6, 1, 10, 3, 6, 1, 6, 3}, // 124 7
    { 0, 1, 9, 10, 7, 6}, // 125 4
    { 7, 0, 8, 7, 6, 0, 3, 0, 10, 10, 0, 6}, // 126 3
    { 7, 6, 10}, // 127 1
    { 7, 10, 6}, // 128 1
    { 3, 8, 0, 10, 6, 7}, // 129 3
    { 0, 9, 1, 10, 6, 7}, // 130 4
    { 8, 9, 1, 8, 1, 3, 10, 6, 7}, // 131 7
    {11, 2, 1, 6, 7, 10}, // 132 3
    { 1, 11, 2, 3, 8, 0, 6, 7, 10}, // 133 6
    { 2, 0, 9, 2, 9, 11, 6, 7, 10}, // 134 7
    { 6, 7, 10, 2, 3, 11, 11, 3, 8, 11, 8, 9}, // 135 12
    { 7, 3, 2, 6, 7, 2}, // 136 2
    { 7, 8, 0, 7, 0, 6, 6, 0, 2}, // 137 5
    { 2, 6, 7, 2, 7, 3, 0, 9, 1}, // 138 7
    { 1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7}, // 139 14
    {11, 6, 7, 11, 7, 1, 1, 7, 3}, // 140 5
    {11, 6, 7, 1, 11, 7, 1, 7, 8, 1, 8, 0}, // 141 9
    { 0, 7, 3, 0, 11, 7, 0, 9, 11, 6, 7, 11}, // 142 11
    { 7, 11, 6, 7, 8, 11, 8, 9, 11}, // 143 5
    { 6, 4, 8, 10, 6, 8}, // 144 2
    { 3, 10, 6, 3, 6, 0, 0, 6, 4}, // 145 5
    { 8, 10, 6, 8, 6, 4, 9, 1, 0}, // 146 7
    { 9, 6, 4, 9, 3, 6, 9, 1, 3, 10, 6, 3}, // 147 11
    { 6, 4, 8, 6, 8, 10, 2, 1, 11}, // 148 7
    { 1, 11, 2, 3, 10, 0, 0, 10, 6, 0, 6, 4}, // 149 12
    { 4, 8, 10, 4, 10, 6, 0, 9, 2, 2, 9, 11}, // 150 10
    {11, 3, 9, 11, 2, 3, 9, 3, 4, 10, 6, 3, 4, 3, 6}, // 151 7
    { 8, 3, 2, 8, 2, 4, 4, 2, 6}, // 152 5
    { 0, 2, 4, 4, 2, 6}, // 153 8
    { 1, 0, 9, 2, 4, 3, 2, 6, 4, 4, 8, 3}, // 154 12
    { 1, 4, 9, 1, 2, 4, 2, 6, 4}, // 155 5
    { 8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 11}, // 156 14
    {11, 0, 1, 11, 6, 0, 6, 4, 0}, // 157 5
    { 4, 3, 6, 4, 8, 3, 6, 3, 11, 0, 9, 3, 11, 3, 9}, // 158 7
    {11, 4, 9, 6, 4, 11}, // 159 2
    { 4, 5, 9, 7, 10, 6}, // 160 3
    { 0, 3, 8, 4, 5, 9, 10, 6, 7}, // 161 6
    { 5, 1, 0, 5, 0, 4, 7, 10, 6}, // 162 7
    {10, 6, 7, 8, 4, 3, 3, 4, 5, 3, 5, 1}, // 163 12
    { 9, 4, 5, 11, 2, 1, 7, 10, 6}, // 164 6
    { 6, 7, 10, 1, 11, 2, 0, 3, 8, 4, 5, 9}, // 165 13
    { 7, 10, 6, 5, 11, 4, 4, 11, 2, 4, 2, 0}, // 166 12
    { 3, 8, 4, 3, 4, 5, 3, 5, 2, 11, 2, 5, 10, 6, 7}, // 167 6
    { 7, 3, 2, 7, 2, 6, 5, 9, 4}, // 168 7
    { 9, 4, 5, 0, 6, 8, 0, 2, 6, 6, 7, 8}, // 169 12
    { 3, 2, 6, 3, 6, 7, 1, 0, 5, 5, 0, 4}, // 170 10
    { 6, 8, 2, 6, 7, 8, 2, 8, 1, 4, 5, 8, 1, 8, 5}, // 171 7
    { 9, 4, 5, 11, 6, 1, 1, 6, 7, 1, 7, 3}, // 172 12
    { 1, 11, 6, 1, 6, 7, 1, 7, 0, 8, 0, 7, 9, 4, 5}, // 173 6
    { 4, 11, 0, 4, 5, 11, 0, 11, 3, 6, 7, 11, 3, 11, 7}, // 174 7
    { 7, 11, 6, 7, 8, 11, 5, 11, 4, 4, 11, 8}, // 175 3
    { 6, 5, 9, 6, 9, 10, 10, 9, 8}, // 176 5
    { 3, 10, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9}, // 177 9
    { 0, 8, 10, 0, 10, 5, 0, 5, 1, 5, 10, 6}, // 178 14
    { 6, 3, 10, 6, 5, 3, 5, 1, 3}, // 179 5
    { 1, 11, 2, 9, 10, 5, 9, 8, 10, 10, 6, 5}, // 180 12
    { 0, 3, 10, 0, 10, 6, 0, 6, 9, 5, 9, 6, 1, 11, 2}, // 181 6
    {10, 5, 8, 10, 6, 5, 8, 5, 0, 11, 2, 5, 0, 5, 2}, // 182 7
    { 6, 3, 10, 6, 5, 3, 2, 3, 11, 11, 3, 5}, // 183 3
    { 5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8}, // 184 11
    { 9, 6, 5, 9, 0, 6, 0, 2, 6}, // 185 5
    { 1, 8, 5, 1, 0, 8, 5, 8, 6, 3, 2, 8, 6, 8, 2}, // 186 7
    { 1, 6, 5, 2, 6, 1}, // 187 2
    { 1, 6, 3, 1, 11, 6, 3, 6, 8, 5, 9, 6, 8, 6, 9}, // 188 7
    {11, 0, 1, 11, 6, 0, 9, 0, 5, 5, 0, 6}, // 189 3
    { 0, 8, 3, 5, 11, 6}, // 190 4
    {11, 6, 5}, // 191 1
    {10, 11, 5, 7, 10, 5}, // 192 2
    {10, 11, 5, 10, 5, 7, 8, 0, 3}, // 193 7
    { 5, 7, 10, 5, 10, 11, 1, 0, 9}, // 194 7
    {11, 5, 7, 11, 7, 10, 9, 1, 8, 8, 1, 3}, // 195 10
    {10, 2, 1, 10, 1, 7, 7, 1, 5}, // 196 5
    { 0, 3, 8, 1, 7, 2, 1, 5, 7, 7, 10, 2}, // 197 12
    { 9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 10}, // 198 14
    { 7, 2, 5, 7, 10, 2, 5, 2, 9, 3, 8, 2, 9, 2, 8}, // 199 7
    { 2, 11, 5, 2, 5, 3, 3, 5, 7}, // 200 5
    { 8, 0, 2, 8, 2, 5, 8, 5, 7, 11, 5, 2}, // 201 11
    { 9, 1, 0, 5, 3, 11, 5, 7, 3, 3, 2, 11}, // 202 12
    { 9, 2, 8, 9, 1, 2, 8, 2, 7, 11, 5, 2, 7, 2, 5}, // 203 7
    { 1, 5, 3, 3, 5, 7}, // 204 8
    { 0, 7, 8, 0, 1, 7, 1, 5, 7}, // 205 5
    { 9, 3, 0, 9, 5, 3, 5, 7, 3}, // 206 5
    { 9, 7, 8, 5, 7, 9}, // 207 2
    { 5, 4, 8, 5, 8, 11, 11, 8, 10}, // 208 5
    { 5, 4, 0, 5, 0, 10, 5, 10, 11, 10, 0, 3}, // 209 14
    { 0, 9, 1, 8, 11, 4, 8, 10, 11, 11, 5, 4}, // 210 12
    {11, 4, 10, 11, 5, 4, 10, 4, 3, 9, 1, 4, 3, 4, 1}, // 211 7
    { 2, 1, 5, 2, 5, 8, 2, 8, 10, 4, 8, 5}, // 212 11
    { 0, 10, 4, 0, 3, 10, 4, 10, 5, 2, 1, 10, 5, 10, 1}, // 213 7
    { 0, 5, 2, 0, 9, 5, 2, 5, 10, 4, 8, 5, 10, 5, 8}, // 214 7
    { 9, 5, 4, 2, 3, 10}, // 215 4
    { 2, 11, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8}, // 216 9
    { 5, 2, 11, 5, 4, 2, 4, 0, 2}, // 217 5
    { 3, 2, 11, 3, 11, 5, 3, 5, 8, 4, 8, 5, 0, 9, 1}, // 218 6
    { 5, 2, 11, 5, 4, 2, 1, 2, 9, 9, 2, 4}, // 219 3
    { 8, 5, 4, 8, 3, 5, 3, 1, 5}, // 220 5
    { 0, 5, 4, 1, 5, 0}, // 221 2
    { 8, 5, 4, 8, 3, 5, 9, 5, 0, 0, 5, 3}, // 222 3
    { 9, 5, 4}, // 223 1
    { 4, 7, 10, 4, 10, 9, 9, 10, 11}, // 224 5
    { 0, 3, 8, 4, 7, 9, 9, 7, 10, 9, 10, 11}, // 225 12
    { 1, 10, 11, 1, 4, 10, 1, 0, 4, 7, 10, 4}, // 226 11
    { 3, 4, 1, 3, 8, 4, 1, 4, 11, 7, 10, 4, 11, 4, 10}, // 227 7
    { 4, 7, 10, 9, 4, 10, 9, 10, 2, 9, 2, 1}, // 228 9
    { 9, 4, 7, 9, 7, 10, 9, 10, 1, 2, 1, 10, 0, 3, 8}, // 229 6
    {10, 4, 7, 10, 2, 4, 2, 0, 4}, // 230 5
    {10, 4, 7, 10, 2, 4, 8, 4, 3, 3, 4, 2}, // 231 3
    { 2, 11, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4}, // 232 14
    { 9, 7, 11, 9, 4, 7, 11, 7, 2, 8, 0, 7, 2, 7, 0}, // 233 7
    { 3, 11, 7, 3, 2, 11, 7, 11, 4, 1, 0, 11, 4, 11, 0}, // 234 7
    { 1, 2, 11, 8, 4, 7}, // 235 4
    { 4, 1, 9, 4, 7, 1, 7, 3, 1}, // 236 5
    { 4, 1, 9, 4, 7, 1, 0, 1, 8, 8, 1, 7}, // 237 3
    { 4, 3, 0, 7, 3, 4}, // 238 2
    { 4, 7, 8}, // 239 1
    { 9, 8, 11, 11, 8, 10}, // 240 8
    { 3, 9, 0, 3, 10, 9, 10, 11, 9}, // 241 5
    { 0, 11, 1, 0, 8, 11, 8, 10, 11}, // 242 5
    { 3, 11, 1, 10, 11, 3}, // 243 2
    { 1, 10, 2, 1, 9, 10, 9, 8, 10}, // 244 5
    { 3, 9, 0, 3, 10, 9, 1, 9, 2, 2, 9, 10}, // 245 3
    { 0, 10, 2, 8, 10, 0}, // 246 2
    { 3, 10, 2}, // 247 1
    { 2, 8, 3, 2, 11, 8, 11, 9, 8}, // 248 5
    { 9, 2, 11, 0, 2, 9}, // 249 2
    { 2, 8, 3, 2, 11, 8, 0, 8, 1, 1, 8, 11}, // 250 3
    { 1, 2, 11}, // 251 1
    { 1, 8, 3, 9, 8, 1}, // 252 2
    { 0, 1, 9}, // 253 1
    { 0, 8, 3}, // 254 1
    {} // 255 0
  };
}
