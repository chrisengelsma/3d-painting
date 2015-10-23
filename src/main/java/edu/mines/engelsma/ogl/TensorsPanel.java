/****************************************************************************
Copyright 2009, Colorado School of Mines and others.
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

import edu.mines.jtk.dsp.EigenTensors3;
import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.sgl.*;
import static edu.mines.jtk.util.ArrayMath.*;
import static edu.mines.jtk.ogl.Gl.GL_AMBIENT_AND_DIFFUSE;

import java.awt.Color;

/**
 * An axis-aligned panel that draws a 3D representation of tensors fields as
 * ellipsoids.
 * The ellipsoid centers are located along a 2D slice within a 3D array. The
 * tensors panel makes use of <em>precomputed</em> eigenvectors and their
 * corresponding eigenvalues in order to construct each ellipsoid by using
 * an ellipsoid glyph.
 * @see EllipsoidGlyph
 * @author Chris Engelsma and Dave Hale, Colorado School of Mines.
 * @version 2009.08.25
 */
public class TensorsPanel extends AxisAlignedPanel {

  /**
   * Constructs a tensor panel with provided spatial coordinates and
   * eigentensor information.
   * @param xyz spatial coordinates.
   * @param et eigentensors information.
   */
  public TensorsPanel(
    Sampling sx, Sampling sy, Sampling sz, EigenTensors3 et)
  {
    _sx = sx;
    _sy = sy;
    _sz = sz;
    _et = et;
    _emax = findMaxEigenvalue();
    setDefaultStates();
  }

  /**
   * Updates the panel.
   */
  public void update() {
    dirtyDraw();
  }

  /**
   * Sets the maximum size of the ellipsoids.
   * The size refers to the number of samples that the ellipsoids span.
   * @param size the ellipsoid size.
   */
  public void setEllipsoidSize(int size) {
    _ellipsoidSize = size;
  }

  /////////////////////////////////////////////////////////////////////////////
  // protected

  protected void draw(DrawContext dc) {
    AxisAlignedFrame aaf = this.getFrame();
    if (aaf==null)
      return;
    Axis axis = aaf.getAxis();
    drawEllipsoids(axis);
  }

  /////////////////////////////////////////////////////////////////////////////
  // private
  private ImagePanelGroup _ipg = null;
  private Sampling _sx,_sy,_sz;
  private EigenTensors3 _et;
  private EllipsoidGlyph _eg = new EllipsoidGlyph();
  private int _ellipsoidSize = 10;
  private float _emax;

  /**
   * Draws the tensors in the form of ellipsoids within the current axis.
   */
  private void drawEllipsoids(Axis axis) {
    // Get sampling information
    int nx = _sx.getCount();
    int ny = _sy.getCount();
    int nz = _sz.getCount();
    double dx = _sx.getDelta();
    double dy = _sy.getDelta();
    double dz = _sz.getDelta();
    double fx = _sx.getFirst();
    double fy = _sy.getFirst();
    double fz = _sz.getFirst();

    // Get the dimensions of the data
    double xmin = _sx.getFirst();
    double xmax = _sx.getLast();
    double ymin = _sy.getFirst();
    double ymax = _sy.getLast();
    double zmin = _sz.getFirst();
    double zmax = _sz.getLast();

    // Maximum length of eigenvectors u, v and w.
    float dmax = 0.5f*(float)_ellipsoidSize;
    float dxmax = (float)dx*dmax;
    float dymax = (float)dy*dmax;
    float dzmax = (float)dz*dmax;

    // Distance in between ellipsoids (in samples).
    // This places spheres (close-to) tangent of one another.
    int kec = (int)(2.0*dmax);

    // Scaling factor for the eigenvectors.
    float scale = dmax/sqrt(_emax);

    // Smallest eigenvalue permitted.
    float etiny = 0.0001f*_emax;

    // Get the current frame.
    AxisAlignedFrame aaf = this.getFrame();


    if (axis==Axis.X) {
      // Y-Axis.
      int nyc = (int)((ymax-ymin)/(2.0f*dymax));
      double dyc = kec*dy;
      double fyc = 0.5f*((ymax-ymin)-(nyc-1)*dyc);
      int jyc = (int)(fyc/dy);

      // Z-Axis.
      int nzc = (int)((zmax-zmin)/(2.0f*dzmax));
      double dzc = kec*dz;
      double fzc = 0.5f*((zmax-zmin)-(nzc-1)*dzc);
      int jzc = (int)(fzc/dz);

      xmin = aaf.getCornerMin().x;
      xmax = aaf.getCornerMax().x;
      ymin = aaf.getCornerMin().y;
      ymax = aaf.getCornerMax().y;
      zmin = aaf.getCornerMin().z;
      zmax = aaf.getCornerMax().z;

      // X-Axis.
      float xc = 0.5f*(float)(xmax+xmin);
      int ix = _sx.indexOfNearest(xc);

      for (int iy=jyc; iy<ny; iy+=kec) {
        float yc = (float)(fy+iy*dy);
        if (ymin<yc-dymax && yc+dymax<ymax) {
          for (int iz=jzc; iz<nz; iz+=kec) {
            float zc = (float)(fz+iz*dz);
            if (zmin<zc-dzmax && zc+dzmax<zmax) {
              float[] e = _et.getEigenvalues(iz,iy,ix);
              float[] u = _et.getEigenvectorU(iz,iy,ix);
              float[] v = _et.getEigenvectorV(iz,iy,ix);
              float[] w = _et.getEigenvectorW(iz,iy,ix);
              float eu = e[0], ev = e[1], ew = e[2];
              if (eu<=etiny) eu = etiny;
              if (ev<=etiny) ev = etiny;
              if (ew<=etiny) ew = etiny;
              float uz = u[0], uy = u[1], ux = u[2];
              float vz = v[0], vy = v[1], vx = v[2];
              float wz = w[0], wy = w[1], wx = w[2];
              float su = scale*sqrt(eu);
              float sv = scale*sqrt(ev);
              float sw = scale*sqrt(ew);
              ux *= su*dx; uy *= su*dy; uz *= su*dz;
              vx *= sv*dx; vy *= sv*dy; vz *= sv*dz;
              wx *= sw*dx; wy *= sw*dy; wz *= sw*dz;
              _eg.draw(xc,yc,zc,ux,uy,uz,vx,vy,vz,wx,wy,wz);
            } 
          }
        }
      }
    } else if (axis==Axis.Y) {
      // X-Axis.
      int nxc = (int)((xmax-xmin)/(2.0f*dmax));
      double dxc = kec*dx;
      double fxc = 0.5f*((xmax-xmin)-(nxc-1)*dxc);
      int jxc = (int)(fxc/dx);

      // Z-Axis.
      int nzc = (int)((zmax-zmin)/(2.0f*dmax));
      double dzc = kec*dz;
      double fzc = 0.5f*((zmax-zmin)-(nzc-1)*dzc);
      int jzc = (int)(fzc/dz);

      xmin = aaf.getCornerMin().x;
      xmax = aaf.getCornerMax().x;
      ymin = aaf.getCornerMin().y;
      ymax = aaf.getCornerMax().y;
      zmin = aaf.getCornerMin().z;
      zmax = aaf.getCornerMax().z;

      // Y-Axis.
      float yc = 0.5f*(float)(ymax+ymin);
      int iy = _sy.indexOfNearest(yc);

      for (int ix=jxc; ix<nx; ix+=kec) {
        float xc = (float)(fx+ix*dx);
        if (xc-dmax>xmin || xc+dmax<xmax) {
          for (int iz=jzc; iz<nz; iz+=kec) {
            float zc = (float)(fz+iz*dz);
            if (zc-dmax>zmin || zc+dmax<zmax) {
              float[] e = _et.getEigenvalues(iz,iy,ix);
              float[] u = _et.getEigenvectorU(iz,iy,ix);
              float[] v = _et.getEigenvectorV(iz,iy,ix);
              float[] w = _et.getEigenvectorW(iz,iy,ix);
              float eu = e[0], ev = e[1], ew = e[2];
              if (eu==0)
                eu = 0.001f*_emax;
              if (ev==0)
                ev = 0.001f*_emax;
              if (ew==0)
                ew = 0.001f*_emax;
              float uz = u[0], uy = u[1], ux = u[2];
              float vz = v[0], vy = v[1], vx = v[2];
              float wz = w[0], wy = w[1], wx = w[2];
              float su = scale*sqrt(eu);
              float sv = scale*sqrt(ev);
              float sw = scale*sqrt(ew);
              ux *= su*dx; uy *= su*dy; uz *= su*dz;
              vx *= sv*dx; vy *= sv*dy; vz *= sv*dz;
              wx *= sw*dx; wy *= sw*dy; wz *= sw*dz;
              if (xc>xmin && xc<xmax && zc>zmin && zc<zmax)
                _eg.draw(xc,yc,zc,ux,uy,uz,vx,vy,vz,wx,wy,wz);
            } 
          }
        }
      }
    } else if (axis==Axis.Z) {
      // X-Axis.
      int nxc = (int)((xmax-xmin)/(2.0f*dmax));
      double dxc = kec*dx;
      double fxc = 0.5f*((xmax-xmin)-(nxc-1)*dxc);
      int jxc = (int)(fxc/dx);

      // Y-Axis.
      int nyc = (int)((ymax-ymin)/(2.0f*dmax));
      double dyc = kec*dy;
      double fyc = 0.5f*((ymax-ymin)-(nyc-1)*dyc);
      int jyc = (int)(fyc/dy);

      xmin = aaf.getCornerMin().x;
      xmax = aaf.getCornerMax().x;
      ymin = aaf.getCornerMin().y;
      ymax = aaf.getCornerMax().y;
      zmin = aaf.getCornerMin().z;
      zmax = aaf.getCornerMax().z;

      // Z-Axis.
      float zc = 0.5f*(float)(zmax+zmin);
      int iz = _sz.indexOfNearest(zc);

      for (int ix=jxc; ix<nx; ix+=kec) {
        float xc = (float)(fx+ix*dx);
        if (xc-dmax>xmin || xc+dmax<xmax) {
          for (int iy=jyc; iy<ny; iy+=kec) {
            float yc = (float)(fy+iy*dy);
            if (yc-dmax>ymin || yc+dmax<ymax) {
              float[] e = _et.getEigenvalues(iz,iy,ix);
              float[] u = _et.getEigenvectorU(iz,iy,ix);
              float[] v = _et.getEigenvectorV(iz,iy,ix);
              float[] w = _et.getEigenvectorW(iz,iy,ix);
              float eu = e[0], ev = e[1], ew = e[2];
              if (eu==0)
                eu = 0.001f*_emax;
              if (ev==0)
                ev = 0.001f*_emax;
              if (ew==0)
                ew = 0.001f*_emax;
              float uz = u[0], uy = u[1], ux = u[2];
              float vz = v[0], vy = v[1], vx = v[2];
              float wz = w[0], wy = w[1], wx = w[2];
              float su = scale*sqrt(eu);
              float sv = scale*sqrt(ev);
              float sw = scale*sqrt(ew);
              ux *= su*dx; uy *= su*dy; uz *= su*dz;
              vx *= sv*dx; vy *= sv*dy; vz *= sv*dz;
              wx *= sw*dx; wy *= sw*dy; wz *= sw*dz;
              if (xc>xmin && xc<xmax && yc>ymin && yc<ymax)
                _eg.draw(xc,yc,zc,ux,uy,uz,vx,vy,vz,wx,wy,wz);
            } 
          }
        }
      }
    }
  }

  /**
   * Finds the largest eigenvalue to be used for scaling.
   */
  private float findMaxEigenvalue() {
    float max = 0;
    int n1 = _et.getN1();
    int n2 = _et.getN2();
    int n3 = _et.getN3();
    float[] e = new float[3];
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          _et.getEigenvalues(i1,i2,i3,e);
          float val = max(e[0],e[1],e[2]);
          if (val>max)
            max=val;
        }
      }
    }
    return max;
  }

  /**
   * Sets the default state for 
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

  private void setDefaultStates() {
    setStates(defaultStateSet(Color.CYAN));
  }
}
