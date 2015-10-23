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
package edu.mines.engelsma.sgl;

import edu.mines.jtk.sgl.*;
import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.ogl.Gl.GL_AMBIENT_AND_DIFFUSE;
import static edu.mines.jtk.util.ArrayMath.*;

import java.awt.*;
import java.util.ArrayList;

/**
 * A group of ellipsoids that display a list of eigentensors.
 * @author Chris Engelsma, Colorado School of Mines.
 * @version 2010.01.05
 */
public class EigenTensorsGroup extends Group {

  /**
   * Constructs a new eigentensors group with given eigentensors.
   * Note: Even sampling is assumed.
   * @param et the eigentensors.
   */
  public EigenTensorsGroup(EigenTensors3 et) {
    this(new Sampling(et.getN1()),
         new Sampling(et.getN2()),
         new Sampling(et.getN3()),
         et);
  }

  /**
   * Constructs a new eigentensors group with specified samplings.
   * @param s1 sampling in the first dimension.
   * @param s2 sampling in the second dimension.
   * @param s3 sampling in the third dimension.
   * @param et the eigentensors.
   */
  public EigenTensorsGroup(
    Sampling s1, Sampling s2, Sampling s3, EigenTensors3 et)
  {
    _sx = s3;
    _sy = s2;
    _sz = s1;
    _et = et;
    _emax = findMaxEigenvalue();
    _etn = new EigenTensorsNode();
    this.addChild(_etn);
    setDefaultStates();
  }

  /**
   * Pulls a tensor from a given coordinate.
   * @param x the x-coordinate.
   * @param y the y-coordinate.
   * @param z the z-coordinate.
   * @param storing true, if storing; false, otherwise.
   */
  public void pullTensor(double x, double y, double z, boolean storing) {
    int ix = _sx.indexOfNearest(x);
    int iy = _sy.indexOfNearest(y);
    int iz = _sz.indexOfNearest(z);
    float xc = (float)x, yc = (float)y, zc = (float)z;
    float scale = 1.0f/sqrt(_emax);
    float[] u = _et.getEigenvectorU(iz,iy,ix);
    float[] v = _et.getEigenvectorV(iz,iy,ix);
    float[] w = _et.getEigenvectorW(iz,iy,ix);
    float[] e = _et.getEigenvalues(iz,iy,ix);
    float eu = e[0], ev = e[1], ew = e[2];
    if (eu<=_etiny) eu = _etiny;
    if (ev<=_etiny) ev = _etiny;
    if (ew<=_etiny) ew = _etiny;
    float uz = u[0], uy = u[1], ux = u[2];
    float vz = v[0], vy = v[1], vx = v[2];
    float wz = w[0], wy = w[1], wx = w[2];
    float su = scale*sqrt(eu);
    float sv = scale*sqrt(ev);
    float sw = scale*sqrt(ew);
    ux *= su; uy *= su; uz *= su;
    vx *= sv; vy *= sv; vz *= sv;
    wx *= sw; wy *= sw; wz *= sw;
    if (storing) elist.add(new Ellipsoid(xc,yc,zc,ux,uy,uz,vx,vy,vz,wx,wy,wz));
    else eTemp = new Ellipsoid(xc,yc,zc,ux,uy,uz,vx,vy,vz,wx,wy,wz);
    _etn.setLists(elist,eTemp);
  }

  /**
   * Clears the temporary ellipsoid (better method?)
   */
  public void clearTempTensor() {
    eTemp = null;
    _etn.setLists(elist,eTemp);
  }

  /**
   * Clears all stored ellipsoids.
   */
  public void clearAll() {
    eTemp = null;
    elist.removeAll(elist);
    _etn.setLists(elist,eTemp);
  }

  /**
   * Sets the size of the ellipsoids.
   * @param size the size of the ellipsoids (in samples).
   */
  public void setEllipsoidSize(float size) {
    _eSize = size;
    update();
  }

  /**
   * Gets the size of the ellipsoids.
   * @return the ellipsoid size.
   */
  public float getEllipsoidSize() {
    return _eSize;
  }

  /**
   * Clears the stored tensor closest to the given coordinates.
   * @param x the x-coordinate.
   * @param y the y-coordinate.
   * @param z the z-coordinate.
   */
  public void clearClosestTensor(double x, double y, double z) {
    float closestDistance = 1e6f;
    if (elist!=null) {
      float distance = closestDistance;
      Ellipsoid closestEllipsoid = null;
      for (Ellipsoid e:elist) {
        distance = (float)sqrt(pow(x-e.xc,2)+pow(y-e.yc,2)+pow(z-e.zc,2));
        if (distance<closestDistance) {
          closestDistance = distance;
          closestEllipsoid = e;
        }
      }
      if (closestDistance<_closestPossible && closestEllipsoid!=null) {
        elist.remove(closestEllipsoid);
      }
    }
    update();
  }

  /**
   * Returns this eigen tensor node.
   * @return this node.
   */
  public Node getNode() {
    return _etn;
  }

  /**
   * Updates the scene graph.
   */
  public void update() {
    dirtyDraw();
  }

  ////////////////////////////////////////////////////////////////////////////
  // private

  /**
   * An ellipsoid.
   * Ellipsoids consist of twelve numbers for proper definition:
   * center coordinates (x,y,z)
   * vectors representing the three major axes directions (x,y,z)x3
   */
  private class Ellipsoid {
    float xc,yc,zc;
    float ux,uy,uz;
    float vx,vy,vz;
    float wx,wy,wz;
    public Ellipsoid(
      float xc, float yc, float zc,
      float ux, float uy, float uz,
      float vx, float vy, float vz,
      float wx, float wy, float wz)
    {
      this.xc = xc; this.yc = yc; this.zc = zc;
      this.ux = ux; this.uy = uy; this.uz = uz;
      this.vx = vx; this.vy = vy; this.vz = vz;
      this.wx = wx; this.wy = wy; this.wz = wz;
    }
  }

  private EigenTensors3 _et; // Stored eigentensors
  private Sampling _sx,_sy,_sz; // Eigentensor sampling.
  private float _emax; // Maximum eigenvalue.
  private float _etiny; // Smallest possible eigenvalue.
  private float _eSize = 5.0f; // Maximum ellipsoid radius.
  private float _closestPossible = _eSize; // Closest allowable ellipsoid

  private EigenTensorsNode _etn;
  private ArrayList<Ellipsoid> elist = new ArrayList<Ellipsoid>();
  private Ellipsoid eTemp;

  private class EigenTensorsNode extends Node {

    /**
     * Constructs a new eigen tensor node with null lists.
     */
    public EigenTensorsNode() {
    }

    /**
     * Sets the values of the lists to be drawn.
     * @param elist the list of stuck ellipsoids.
     * @param eTemp the current temporary ellipsoid.
     */
    public void setLists(ArrayList<Ellipsoid> elist, Ellipsoid eTemp) {
      this.elist = elist;
      this.eTemp = eTemp;
      dirtyDraw();
    }

    /**
     * Redraws the scene graph.
     * @param dc the draw context.
     */
    protected void draw(DrawContext dc) {
      if (elist!=null) {
        for (Ellipsoid e:elist)
          _eg.draw(e.xc,e.yc,e.zc,
                   _eSize*e.ux,_eSize*e.uy,_eSize*e.uz,
                   _eSize*e.vx,_eSize*e.vy,_eSize*e.vz,
                   _eSize*e.wx,_eSize*e.wy,_eSize*e.wz);
      }
      if (eTemp!=null)
        _eg.draw(eTemp.xc,eTemp.yc,eTemp.zc,
                 _eSize*eTemp.ux,_eSize*eTemp.uy,_eSize*eTemp.uz,
                 _eSize*eTemp.vx,_eSize*eTemp.vy,_eSize*eTemp.vz,
                 _eSize*eTemp.wx,_eSize*eTemp.wy,_eSize*eTemp.wz);
    }

    private ArrayList<Ellipsoid> elist; // The running list of ellipsoids
    private Ellipsoid eTemp; // A temporary ellipsoid.
    private EllipsoidGlyph _eg = new EllipsoidGlyph();
  }

  /**
   * Sets the state of the group.
   * @param color a color.
   * @return the state set.
   */
  private static StateSet defaultStateSet(Color color) {
    StateSet states = new StateSet();
    ColorState cs = new ColorState();
    cs.setColor(color);
    LightModelState lms = new LightModelState();
    lms.setTwoSide(true);
    MaterialState ms = new MaterialState();
    ms.setColorMaterial(GL_AMBIENT_AND_DIFFUSE);
    ms.setSpecular(Color.WHITE);
    ms.setShininess(100.0f);
    states.add(cs);
    states.add(lms);
    states.add(ms);
    return states;
  }

  /**
   * Sets the default state of the group.
   */
  private void setDefaultStates() {
    setStates(defaultStateSet(Color.CYAN));
  }

  /**
   * Determines the maximum eigenvalue for proper scaling.
   * @return the maximum eigenvalue.
   */
  private float findMaxEigenvalue() {
    int n1 = _et.getN1();
    int n2 = _et.getN2();
    int n3 = _et.getN3();
    float[] e = new float[3];
    float emax = 0.0f;
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          _et.getEigenvalues(i1,i2,i3,e);
          float emaxi = max(e[0],e[1],e[2]);
          if (emax<emaxi)
            emax = emaxi;
        }
      }
    }
    _etiny = 0.0001f*emax;
    return emax;
  }
}
