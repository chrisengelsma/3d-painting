/**
 * Vertex shader for depth-dependent halos
 *
 * Author: Chris Engelsma, Colorado School of Mines
 * Version: 2010.03.23
 */

varying vec2 coords;

void main() {
  coords = gl_Vertex.xy;
  gl_Position = ftransform();
}
