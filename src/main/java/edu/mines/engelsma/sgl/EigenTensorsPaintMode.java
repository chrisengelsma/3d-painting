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
import edu.mines.jtk.awt.Mode;
import edu.mines.jtk.awt.ModeManager;

import java.awt.Component;
import java.awt.event.*;
import java.util.Iterator;
import javax.swing.*;

/**
 * A mode for painting ellipsoids representing eigen-tensors.
 * @author Chris Engelsma, Colorado School of Mines.
 * @version 2010.01.05
 */
public class EigenTensorsPaintMode extends Mode {
  private static final long serialVersionUID = 1L;

  /**
   * Constructs an eigen-tensors edu.mines.engelsma.paint mode with specified 
   * manager.
   * @param modeManager the mode manager for this mode.
   */
  public EigenTensorsPaintMode(ModeManager modeManager) {
    super(modeManager);
    setName("Paint");
    Class<EigenTensorsPaintMode> cls = EigenTensorsPaintMode.class;
    setIcon(loadIcon(cls,"resources/EigenTensorsPaintIcon16.png"));
    setCursor(loadCursor(cls,"resources/EigenTensorsPaintCursor16.png",1,1));
    setMnemonicKey(KeyEvent.VK_P);
    setAcceleratorKey(KeyStroke.getKeyStroke(KeyEvent.VK_P,0));
    setShortDescription("Paint Tensors");
  }

  ///////////////////////////////////////////////////////////////////////////
  // protected

  protected void setActive(Component component, boolean active) {
    if (component instanceof ViewCanvas) {
      if (active) {
        component.addMouseListener(_ml);
        component.addKeyListener(_kl);
        component.addMouseWheelListener(_mwl);
      } else {
        component.removeMouseListener(_ml);
        component.removeKeyListener(_kl);
        component.removeMouseWheelListener(_mwl);
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private ViewCanvas _canvas; // canvas when mouse pressed
  private View _view; // view when mouse pressed; null, if none
  private World _world; // world when mouse pressed; null, if none
  private PickResult _pickResult; // pick result when mouse is pressed.
  private EigenTensorsGroup _etg; // eigen-tensors
  private boolean _sticking; // true used iff sticking ellipsoids.
  private boolean _erasing; // true iff erasing ellipsoids.

  private KeyListener _kl = new KeyListener() {
    public void keyPressed(KeyEvent e) {
      if (e.getKeyChar()=='c') {
        _etg.clearAll();
      }
    }

    // Not used
    public void keyTyped(KeyEvent e) {}
    public void keyReleased(KeyEvent e) {}
  };

  private MouseListener _ml = new MouseAdapter() {

    public void mousePressed(MouseEvent e) {

      // We first assume that we are not sticking or erasing ellipsoids
      _sticking = false;
      _erasing = false;

      // Pick and look in the result for the eigen-tensors group
      _pickResult = pick(e,true);
      _canvas = (ViewCanvas)e.getSource();
      _view = _canvas.getView();
      if (_view!=null)
        _world = _view.getWorld();

      // The check for eigen-tensors group...
      if (_pickResult!=null) {
        Class<EigenTensorsGroup> clsetg = EigenTensorsGroup.class;
        Iterator<Node> itn = _world.getChildren();
        while(itn.hasNext()) {
          Node n = itn.next();
          if (n.getClass()==clsetg) {
            _etg = (EigenTensorsGroup)n;
          }
        }

        if (_etg!=null) {
          Point3 p = _pickResult.getPointLocal();

          // If shift selected, stick the ellipsoids.
          if (e.isShiftDown()) _sticking = true;
          else _sticking = false;

          // If control selected, erase ellipsoids.
          if (e.isAltDown()) _erasing = true;
          else _erasing = false;

          if (_erasing) _etg.clearClosestTensor(p.x,p.y,p.z);
          else _etg.pullTensor(p.x,p.y,p.z,_sticking);
        }
      }
      // Begin listening for mouse movement.
      _canvas.addMouseMotionListener(_mml);
    }

    public void mouseReleased(MouseEvent e) {
      // Make sure the temporary ellipsoid returns to null
      _etg.clearTempTensor();

      // No longer painting.
      _canvas.removeMouseMotionListener(_mml);
    }
  };

  // The mouse wheel changes the size of all the ellipsoids.

  private MouseWheelListener _mwl = new MouseWheelListener() {
    public void mouseWheelMoved(MouseWheelEvent e) {
      // Measure the wheel rotations and scale by +/- rotations^-1.
      if (_etg!=null) {
        int nclicks = e.getWheelRotation();
        float size = _etg.getEllipsoidSize();
        if (!(size<0.5f && nclicks<0)) {
          float scroll = (float)nclicks/10.0f;
          size += scroll;
          _etg.setEllipsoidSize(size);
        }
      }
    }
  };

  private MouseMotionListener _mml = new MouseMotionAdapter() {
    public void mouseDragged(MouseEvent e) {
      if (_etg!=null) {
        // Same rules apply as mousePressed
        _pickResult = pick(e,false);
        if (_pickResult!=null) {
          Point3 p = _pickResult.getPointLocal();
          if (e.isShiftDown()) _sticking = true;
          else _sticking = false;

          if (e.isAltDown()) _erasing = true;
          else _erasing = false;

          if (_erasing) _etg.clearClosestTensor(p.x,p.y,p.z);
          else _etg.pullTensor(p.x,p.y,p.z,_sticking);
        }
      }
    }
  };

  private PickResult pick(MouseEvent event, boolean loud) {
    ViewCanvas canvas = (ViewCanvas)event.getSource();
    View view = canvas.getView();
    if (view==null)
      return null;
    World world = view.getWorld();
    if (world==null)
      return null;
    PickContext pc = new PickContext(event);
//    world.pickApply(pc);
    PickResult pickResult = pc.getClosest();
    if (pickResult!=null) {
      Point3 pointLocal = pickResult.getPointLocal();
      Point3 pointWorld = pickResult.getPointWorld();
      if (loud) {
        System.out.println("Paint Pick");
        System.out.println("  local="+pointLocal);
        System.out.println("  world="+pointWorld);
      }
    } else {
      if (loud)
        System.out.println("Paint Pick nothing");
    }
    return pickResult;
  }
}
