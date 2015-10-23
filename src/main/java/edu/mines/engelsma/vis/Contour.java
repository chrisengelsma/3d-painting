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
package edu.mines.engelsma.vis;

/**
 * A contour.
 * Contours are represented by lists of vertices and triangles.
 * @author Chris Engelsma, Colorado School of Mines
 * @version 2010.12.10
 */
public class Contour {
  /**
   * Vertex array with packed coordinates (x1,x2,x3).
   * The number of vertices equals x.length/3.
   */
  public float[] x;

  /**
   * Normal array with packed components (u1,u2,u3).
   * The number of normal vectors equals u.length/3. This number equals the
   * number of vertices.
   */
  public float[] u;

  /**
   * Triangle array of packed vertex indices (i1,i2,i3).
   * When multiplied by 3, each index references the first coordinate x1 of
   * a vertex (x1,x2,x3) stored in the packed vertex array. The number of
   * triangles equals i.length/3. A vertex may be referenced by more than
   * one triangle.
   */
  public int[] i;
}
