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
package edu.mines.engelsma.ogl;

import edu.mines.jtk.bench.ArrayListBench.FloatList;
import static edu.mines.jtk.util.ArrayMath.*;
import edu.mines.jtk.dsp.EigenTensors3;

/**
 * Draws an ellipsoid given major axes values and associated scaling factors.
 * @author Chris Engelsma, Colorado School of Mines.
 * @version 2010.02.23
 */
public class Ellipsoid {

  /**
   * Constructs an ellipsoid given three orthonormal vectors (eigenvectors) 
   * and scaling factors (eigenvalues).
   * @param u a normalized eigenvector (ux,uy,uz).
   * @param v a normalized eigenvector (vx,vy,vz).
   * @param w a normalized eigenvector (wx,wy,wz).
   * @param a a float array of three eigenvalues.
   */
  public Ellipsoid(float[] u, float[] v, float[] w, float[] a) {
    if (v!=null) {
      _v = v;
    } else {
      _v = normCrossProd(u,w);
    }
    _u = u;
    _w = w;
    _a = a;
    _A = getEigenTensor();
    calculateVertices();
  }

  /**
   * Constructs an ellipsoid given two orthonormal vectors (eigenvectors) and
   * three scaling factors (eigenvalues).
   * The third orthonormal vector is computed by using the normalized cross
   * product of the two given vectors.
   * @param u a normalized eigenvector (ux,uy,uz).
   * @param w a normalized eigenvector (vx,vy,vz).
   * @param a a float array of three eigenvalues.
   */
  public Ellipsoid(float[] u, float[] w, float[] a) {
    this(u,null,w,a);
  }

  /**
   * Gets the triangle vertices used to render the ellipsoid.
   * @return the triangle vertices.
   */
  public float[] getVertices() {
    return _vertices.trim();
  }

  /**
   * Gets the eigenvector U(x,y,z).
   * @return the eigenvector U.
   */
  public float[] getEigenvectorU() {
    return _u;
  }

  /**
   * Gets the eigenvector V(x,y,z).
   * @return the eigenvector V.
   */
  public float[] getEigenvectorV() {
    return _v;
  }

  /**
   * Gets the eigenvector W(x,y,z).
   * @return the eigenvector W.
   */
  public float[] getEigenvectorW() {
    return _w;
  }

  //////////////////////////////////////////////////////////////////////////////
  // private

  private static float[] _u;
  private static float[] _v;
  private static float[] _w;
  private static float[] _a;
  private static float[][] _A;
  private static FloatList _vertices = new FloatList();

  private float[][] getEigenTensor() {
    int n3,n2,n1;
    n3 = n2 = n1 = 1;
    EigenTensors3 et = new EigenTensors3(n3,n2,n1,false);
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          et.setEigenvalues(i1,i2,i3,_a);
          et.setEigenvectorU(i1,i2,i3,_u);
          et.setEigenvectorW(i1,i2,i3,_w);
        }
      }
    }
    float[] t = et.getTensor(0,0,0);
    float[][] A = new float[3][3];
    A[0][0] = t[0];
    A[1][1] = t[3];
    A[2][2] = t[5];
    A[0][1] = A[1][0] = t[1];
    A[0][2] = A[2][0] = t[2];
    A[1][2] = A[2][1] = t[4];

    printTensor(A);
    return A;
  }

  /**
   * DELETE
   */
  private void printTensor(float[][] t) {
    for (int i=0; i<3; ++i) {
      for (int j=0; j<3; ++j) {
        System.out.print(t[i][j]+"\t");
      }
      System.out.println();
    }
  }

  private float[] normCrossProd(float[] u, float[] w) {
    float a = u[1]*w[2]-u[2]*w[1];
    float b = u[2]*w[0]-u[0]*w[2];
    float c = u[0]*w[1]-u[1]*u[0];
    float s = 1/sqrt(a*a+b*b+c*c);
    return new float[]{a*s,b*s,c*s};
  }

  /**
   * Calculates the vertices for an icosahedron then subdivides.
   */
  private void calculateVertices() {

    /*
    float tau = (1.0f + sqrt(5.0f))/2.0f; // golden ratio

    float X = _orthoVector[0];
    float Y = _orthoVector[1];
    float Z = _orthoVector[2];

    float vX = X/(sqrt(tau*sqrt(5)));
    float vY = Y/(sqrt(tau*sqrt(5)));
    float vZ = Z/(sqrt(tau*sqrt(5)));


    float[][] vertices =
      {{-vX,0.0f,tau*vZ},{vX,0.0f,tau*vZ},    // XZ-Plane
       {-vX,0.0f,-tau*vZ},{vX,0.0f,-tau*vZ},
       {0.0f,tau*vY,vZ},{0.0f,tau*vY,-vZ},    // YZ-Plane
       {0.0f,-tau*vY,vZ},{0.0f,-tau*vY,-vZ},
       {tau*vX,vY,0.0f},{-tau*vX,vY,0.0f},    // XY-Plane
       {tau*vX,-vY,0.0f},{-tau*vX,-vY,0.0f}};

    int[] indices =
      {1,4,0,     4,9,0,    4,5,9,    8,5,4,    1,8,4,
       1,10,8,    10,3,8,   8,3,5,    3,2,5,    3,7,2,
       3,10,7,    10,6,7,   6,11,7,   6,0,11,   6,1,0,
       10,1,6,    11,0,9,   2,11,9,   5,2,9,    11,2,7};

    float[] icov = new float[indices.length*3]; // icosahedron vertices

    for (int i=0; i<indices.length; i++) {
      icov[(3*i)+0] = vertices[indices[i]][0];
      icov[(3*i)+1] = vertices[indices[i]][1];
      icov[(3*i)+2] = vertices[indices[i]][2];
    }
     */

    float[] u = _u;
    float[] v = _v;
    float[] w = _w;

    for (int i=0; i<3; ++i) {
      _u[i] *= _a[0];
      _v[i] *= _a[1];
      _w[i] *= _a[2];
    }

    float[][] vertices = 
      {{u[0],u[1],u[2]},{-u[0],-u[1],-u[2]},
       {v[0],v[1],v[2]},{-v[0],-v[1],-v[2]},
       {w[0],w[1],w[2]},{-w[0],-w[1],-w[2]}};

    int[] tindices =
      {4,0,3,  4,3,1,  4,1,2,  4,2,0,
       5,0,2,  5,2,1,  5,1,3,  5,3,0};

    float[] icov = new float[tindices.length*3];

    for (int i=0; i<tindices.length; i++) {
      icov[(3*i)+0] = vertices[tindices[i]][0];
      icov[(3*i)+1] = vertices[tindices[i]][1];
      icov[(3*i)+2] = vertices[tindices[i]][2];
    }

    for (int i=0; i<icov.length-8; i+=9) {
      float[] v1 = new float[]
        {icov[i+0],icov[i+1],icov[i+2]};
      float[] v2 = new float[]
        {icov[i+3],icov[i+4],icov[i+5]};
      float[] v3 = new float[]
        {icov[i+6],icov[i+7],icov[i+8]};
      subdivide(v1,v2,v3,6);
    }
  }

  /**
   * Subdivides the segments of a given triangle.
   */
  private void subdivide(
    float[] v1, float[] v2, float[] v3, float depth)
  {
    float[] v12 = new float[3];
    float[] v23 = new float[3];
    float[] v31 = new float[3];


    // If done, add the final vertices and return.
    if (depth==0) {
      for (int i=0; i<3; i++)
        _vertices.add(v1[i]);
      for (int i=0; i<3; i++)
        _vertices.add(v2[i]);
      for (int i=0; i<3; i++)
        _vertices.add(v3[i]);
      return;
    }

    // Find the midpoint.
    for (int i=0; i<3; i++) {
      v12[i] = (v1[i]+v2[i])/2.0f;
      v23[i] = (v2[i]+v3[i])/2.0f;
      v31[i] = (v3[i]+v1[i])/2.0f;
    }

    // Normalize the coordinates.
    v12 = normalize(v12);
    v23 = normalize(v23);
    v31 = normalize(v31);

    // Recursively subdivide
    subdivide(v2,v23,v12,depth-1);
    subdivide(v1,v12,v31,depth-1);
    subdivide(v3,v31,v23,depth-1);
    subdivide(v12,v23,v31,depth-1);
  }

  private float[] normalize(float[] v) {
    float x = v[0];
    float y = v[1];
    float z = v[2];
    float d = sqrt(x*x+y*y+z*z);
    v[0] /= d;
    v[1] /= d;
    v[2] /= d;
    return v;
  }

  /*
   * Normalizes the given vector and applies the appropriate
   * radius.
   */
  /*
  private float[] normalize(float[] v) {
    float x = v[0];
    float y = v[1];
    float z = v[2];

    float[] vn = new float[3];

    // Grad[F(x,y,z)]
    vn[0] = 2*(_A[0][0]*x+_A[0][1]*y+_A[0][2]*z);
    vn[1] = 2*(_A[1][0]*x+_A[1][1]*y+_A[1][2]*z);
    vn[2] = 2*(_A[2][0]*x+_A[2][1]*y+_A[2][2]*z);

    // Scale by 1/sqrt(x'Ax)
    float xAx = xAx(x,y,z);
    vn[0] /= sqrt(xAx);
    vn[1] /= sqrt(xAx);
    vn[2] /= sqrt(xAx);

    return vn;
  }
   */

  private float xAx(float x, float y, float z) {
    return (_A[0][0]*x*x+_A[1][1]*y*y+_A[2][2]*z*z +
            2*(_A[0][1]*x*y+_A[0][2]*x*z+_A[1][2]*y*z));
  }
}
