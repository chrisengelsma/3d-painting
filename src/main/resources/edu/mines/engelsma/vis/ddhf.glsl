/**
 * Fragment (pixel) shader for depth-dependent halos
 *
 * Author: Chris Engelsma, Colorado School of Mines
 * Version: 2010.03.23
 */

varying vec2 coords;
uniform sampler2D tex;

void main() {
  gl_FragColor = texture2D(tex,coords);
}
