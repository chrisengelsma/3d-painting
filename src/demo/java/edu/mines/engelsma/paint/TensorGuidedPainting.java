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

import edu.mines.engelsma.util.FakeData;
import edu.mines.engelsma.vis.Contour;

import edu.mines.jtk.awt.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.io.*;
import edu.mines.jtk.sgl.*;

import static edu.mines.jtk.ogl.Gl.*;
import static edu.mines.jtk.util.ArrayMath.*;

import javax.swing.*;
import javax.swing.event.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.image.*;
import java.io.*;

/**
 * Demos 3D painting.
 * The default demo constructs a 3D function that the user can
 * edu.mines.engelsma.paint. The
 * user also has the option to load in an image and precomputed tensors.
 * @author Chris Engelsma, Colorado School of Mines
 * @version 2010.12.10
 */
public class TensorGuidedPainting {

  /**
   * Constructs a new tensor-guided painting demo.
   * @param n1 size of the 1st dimension.
   * @param n2 size of the 2nd dimension.
   * @param n3 size of the 3rd dimension.
   * @param dm a distance map.
   * @param image a 3D image.
   */
  public TensorGuidedPainting(
    int n1, int n2, int n3,
    PaintBrush dm, float[][][] image)
  {
    this(new Sampling(n1),new Sampling(n2),new Sampling(n3),dm,image);
  }

  /**
   * Constructs a new tensor-guided painting demo.
   * @param s1 sampling in the 1st dimension.
   * @param s2 sampling in the 2nd dimension.
   * @param s3 sampling in the 3rd dimension.
   * @param dm a distance map.
   * @param image a 3D image.
   */
  public TensorGuidedPainting(
    Sampling s1, Sampling s2, Sampling s3,
    PaintBrush dm, float[][][] image)
  {
    _s1 = s1;
    _s2 = s2;
    _s3 = s3;
    _n1 = s1.getCount();
    _n2 = s2.getCount();
    _n3 = s3.getCount();
    _dm = dm;
    _image = image;
    
    float[][][] paint = new float[_n3][_n2][_n1];
    _p3g = new Painting3Group(s1,s2,s3,paint);
    _p3 = _p3g.getPainting3();
    _ipg = new ImagePanelGroup2(s1,s2,s3,image,_p3g.getPaint());
    _ipg.setColorModel1(ColorMap.getGray());
    _ipg.setColorModel2(_icm);
    _ipg.setClips2(0.0f,1.0f);

    _frame = new SimpleFrame();
    _dm.setSize(30);
    _cm = new ColorMap(_icm);
    _pc = new PaintControl(_frame);
    
    ModeManager mm = _frame.getModeManager();
    PaintBrushMode pbm = new PaintBrushMode(mm);

    JToolBar jtb = _frame.getJToolBar();
    JToggleButton pbmButton = new ModeToggleButton(pbm);
    jtb.add(pbmButton);

    _world = _frame.getWorld();
    _world.addChild(_ipg);
    _world.addChild(_p3g);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private class PaintControl extends JDialog {
    public PaintControl(JFrame parent) {
      super(parent,false);
      setAlwaysOnTop(true);
      JPanel panel = new JPanel();
      panel.setLayout(new BoxLayout(panel,BoxLayout.Y_AXIS));

      _ac = new AlphaChooser(this);
      setAlpha(_alpha);
      JButton btn = new JButton("Toggle rendering");
      btn.addActionListener(new ActionListener() {
        public void actionPerformed(ActionEvent e) {
          if (rendered) {
            remove3DRendering();
            rendered = false;
          } else {
            renderPaintedVoxels(_brushColor,amplitudes);
            rendered = true;
          }
        }
      });

      JButton btn2 = new JButton("Toggle amplitudes on volume");
      btn2.addActionListener(new ActionListener() {
        public void actionPerformed(ActionEvent e) {
          if (amplitudes) amplitudes = false;
          else amplitudes = true;
          if (rendered) renderPaintedVoxels(_brushColor,amplitudes);
        }
      });

      panel.add(new JLabel("alpha"));
      panel.add(_ac);
      panel.add(btn);
      panel.add(btn2);

      this.add(panel);
      this.pack();
      this.setVisible(true);
      this.setLocation(new Point(_frame.getWidth(),0));
    }
    void setAlpha(double alpha) {
      _alpha = alpha;
      updatePaintColorModel(alpha);
    }
    double getAlpha() {
      return _alpha;
    }
    private AlphaChooser _ac;
    private double _alpha = 0.3;
    private boolean rendered = false;
    private boolean amplitudes = false;
  }

  private class AlphaChooser extends JPanel {
    public AlphaChooser(TensorGuidedPainting.PaintControl pc) {
      super();
      _pc = pc;
      _slider = new JSlider(0,100);
      _slider.setValue((int)(100*_pc.getAlpha()));
      _slider.addChangeListener(new ChangeListener() {
        public void stateChanged(ChangeEvent e) {
          if (_slider==e.getSource()) {
            int value = _slider.getValue();
            _pc.setAlpha(0.01*(double)value);
          }
      }
      });
      this.add(_slider);
    }
    private JSlider _slider;
    private TensorGuidedPainting.PaintControl _pc;
  }
  
  private class PaintBrushMode extends Mode {
    public PaintBrushMode(ModeManager modeManager) {
      super(modeManager);
      setName("Paint");
      Class<PaintBrushMode> cls = PaintBrushMode.class;
      setIcon(loadIcon(cls,"PaintIcon16.png"));
      setCursor(loadCursor(cls,"PaintCursor16.png",1,1));
      setMnemonicKey(KeyEvent.VK_B);
      setAcceleratorKey(KeyStroke.getKeyStroke(KeyEvent.VK_P,0));
      setShortDescription("Paint brush");
    }

    /////////////////////////////////////////////////////////////////////////
    // protected

    protected void setActive(Component component, boolean active) {
      if (component instanceof ViewCanvas) {
        if (active) {
          component.addMouseListener(_ml);
          component.addMouseWheelListener(_mwl);
          component.addKeyListener(_kl);
        } else {
          component.removeMouseListener(_ml);
          component.removeMouseWheelListener(_mwl);
          component.removeKeyListener(_kl);
        }
      }
    }

    /////////////////////////////////////////////////////////////////////////
    // private

    private ViewCanvas _canvas;
    private World _world;
    private boolean _painting;
    private boolean _erasing = false;
    private MouseConstrained _mouseConstrained;
    private TriangleGroup _dmtg = null;

    private KeyListener _kl = new KeyListener() {
      public void keyTyped(KeyEvent e) {}
      public void keyPressed(KeyEvent e) {}
      public void keyReleased(KeyEvent e) {
        if (e.getKeyChar()=='e') {
          if (_erasing) _erasing=false;
          else _erasing = true;
          updateContour();
          updatePaint();
        }
      }
    };

    private MouseListener _ml = new MouseAdapter() {
      public void mousePressed(MouseEvent e) {
        PickResult pr = pick(e);
        if (pr!=null) {
          Node node = pr.getNode(AxisAlignedQuad.class);
          if (node!=null) {
            _painting = true;
            _canvas.addMouseMotionListener(_mml);
            AxisAlignedQuad quad = (AxisAlignedQuad)node;
            AxisAlignedFrame frame = quad.getFrame();
            Axis axis = frame.getAxis();
            Point3 origin = pr.getPointWorld();
            Vector3 normal = null;
            if (axis==Axis.X) {
              normal = new Vector3(1.0,0.0,0.0);
            } else if (axis==Axis.Y) {
              normal = new Vector3(0.0,1.0,0.0);
            } else if (axis==Axis.Z) {
              normal = new Vector3(0.0,0.0,1.0);
            }
            Plane plane = new Plane(origin,normal);
            Matrix44 worldToPixel = pr.getWorldToPixel();
            _mouseConstrained = new MouseOnPlane(e,origin,plane,worldToPixel);
            paintAt(origin);
          }
        }
      }
      public void mouseReleased(MouseEvent e) {
        if (_painting) {
          _mouseConstrained = null;
          _canvas.removeMouseMotionListener(_mml);
          _painting = false;
        }
      }
    };
    private MouseMotionListener _mml = new MouseMotionAdapter() {
      public void mouseDragged(MouseEvent e) {
        Point3 point = _mouseConstrained.getPoint(e);
        paintAt(point);
      }
    };
    private MouseWheelListener _mwl = new MouseWheelListener() {
      public void mouseWheelMoved(MouseWheelEvent e) {
        int nclicks = e.getWheelRotation();
        int size = _dm.getSize();
        size += nclicks;
        _dm.setSize(size);
        if (_dmtg==null) return;
        updateContour();
        updatePaint();
      }
    };

    private void paintAt(Point3 point) {
      int i1 = max(0,min(_n1-1,(int)(_s1.indexOfNearest(point.z)+0.5)));
      int i2 = max(0,min(_n2-1,(int)(_s2.indexOfNearest(point.y)+0.5)));
      int i3 = max(0,min(_n3-1,(int)(_s3.indexOfNearest(point.x)+0.5)));
      paintAt(i1,i2,i3);
    }

    private void paintAt(int i1, int i2, int i3) {
      _dm.setLocation(i1,i2,i3);
      updateContour();
      updatePaint(i1,i2,i3);
    }

    private void updateContour() {
      if (_dmtg!=null)
        _world.removeChild(_dmtg);
      Contour contour = _dm.getContour();
      _dmtg = new TriangleGroup(contour.i,contour.x,contour.u);
      StateSet states = new StateSet();
      ColorState cs = new ColorState();
      ColorMap cm = new ColorMap(_ipg.getColorModel2());
      Color color = cm.getColor(_brushColor);
      if (_erasing) cs.setColor(Color.BLACK);
      else cs.setColor(color);
      states.add(cs);
      LightModelState lms = new LightModelState();
      lms.setTwoSide(true);
      states.add(lms);
      MaterialState ms = new MaterialState();
      ms.setColorMaterial(GL_AMBIENT_AND_DIFFUSE);
      ms.setSpecular(Color.WHITE);
      ms.setShininess(100.0f);
      states.add(ms);
      _dmtg.setStates(states);
      _world.addChild(_dmtg);
    }

    private void updatePaint() {
      int[] k = _dm.getLocation();
      updatePaint(k[0],k[1],k[2]);
    }

    private void updatePaint(int i1, int i2, int i3) {
      float size = _dm.getSize();
      if (_erasing) _p3.eraseAt(i1,i2,i3,size,_dm);
      else          _p3.paintAt(i1,i2,i3,_brushColor,size,_dm);
      _ipg.update2();
    }

    private PickResult pick(MouseEvent event) {
      _canvas = (ViewCanvas)event.getSource();
      View view = _canvas.getView();
      if (view==null)
        return null;
      _world = view.getWorld();
      if (_world==null)
        return null;
      PickContext pc = new PickContext(event);
      _world.pick(pc);
      PickResult pickResult = pc.getClosest();
      if (pickResult!=null) {
        Point3 pointLocal = pickResult.getPointLocal();
        Point3 pointWorld = pickResult.getPointWorld();
        System.out.println("Painting at:");
        System.out.println("  local="+pointLocal);
        System.out.println("  world="+pointWorld);
      } else {
        System.out.println("Painting nothing");
      }

      return pickResult;
    }
  }

  private World _world;
  private int _n1,_n2,_n3;
  private Sampling _s1,_s2,_s3;
  private SimpleFrame _frame;
  private ColorMap _cm;
  private ImagePanelGroup2 _ipg;
  private PaintControl _pc;
  private PaintBrush _dm;
  private Painting3 _p3;
  private Painting3Group _p3g;
  private IndexColorModel _icm = ColorMap.getJet();
  private float[][][] _image;
  private TriangleGroup _tg;

  private static EigenTensors3 _et;
  private static float _brushColor = 0.5f;

  private void remove3DRendering() {
    _world.removeChild(_tg);
  }

  private void renderPaintedVoxels(float value, boolean amplitudes) {
    if (_tg!=null) _world.removeChild(_tg);
    float[] rgb;
    Contour c = _p3.getContour(value);
    if (amplitudes) {
      rgb = extractImageColors(c.x);
      _tg = new TriangleGroup(c.i,c.x,null,rgb);
    } else
      _tg = new TriangleGroup(c.i,c.x);
    setModelState(_tg);
    _world.addChild(_tg);
  }

  private void setModelState(TriangleGroup tg) {
    StateSet states = new StateSet();
    LightModelState lms = new LightModelState();
    lms.setTwoSide(true);
    ColorState cs = new ColorState();
    ColorMap cm = new ColorMap(_ipg.getColorModel2());
    Color color = cm.getColor(_brushColor);
    cs.setColor(color);
    MaterialState ms = new MaterialState();
    ms.setColorMaterial(GL_AMBIENT_AND_DIFFUSE);
    ms.setShininess(0.0f);
    states.add(lms);
    states.add(cs);
    states.add(ms);
    tg.setStates(states);
  }

  private float[] extractImageColors(float[] v) {
    int i1j,i2j,i3j;
    int i1k,i2k,i3k;
    double p,v1,v2,v3,val;
    double d1 = _s1.getDelta(), d2 = _s2.getDelta(), d3 = _s3.getDelta();
    Color c;
    ColorMap cm = new ColorMap(_ipg.getColorModel1());
    cm.setValueRange(_ipg.getClip1Min(),_ipg.getClip1Max());
    float[] vals = new float[v.length];
    for (int i=0; i<v.length-2; i+=3) {
      v3 = v[i+0];
      v2 = v[i+1];
      v1 = v[i+2];
      i3j = _s3.indexOfNearest(v3-d3/2);
      i2j = _s2.indexOfNearest(v2-d2/2);
      i1j = _s1.indexOfNearest(v1-d1/2);
      byte[] ins = _p3.getEdgeIntersectionsAt(i1j,i2j,i3j);
      for (int j=0; j<3; ++j) {
        if (ins[j]!=-1) {
          switch(j) {
            case(0):
              i1k = _s1.indexOfNearest(v1+d1/2);
              i2k = i2j;
              i3k = i3j;
              break;
            case(1):
              i1k = i1j;
              i2k = _s2.indexOfNearest(v2+d2/2);
              i3k = i3j;
              break;
            default:
              i1k = i1j;
              i2k = i2j;
              i3k = _s3.indexOfNearest(v3+d3/2);
          }
          p = (double)ins[j]/100.0;
          val = (1.0-p)*_image[i3j][i2j][i1j]+p*_image[i3k][i2k][i1k];
          c = cm.getColor(val);
          vals[i+0] = c.getRed()/255.0f;
          vals[i+1] = c.getGreen()/255.0f;
          vals[i+2] = c.getBlue()/255.0f;
        } else if (ins[j]==-1 && j==2) {
          val = _image[i3j][i2j][i1j];
          c = cm.getColor(val);
          vals[i+0] = c.getRed()/255.0f;
          vals[i+1] = c.getRed()/255.0f;
          vals[i+2] = c.getRed()/255.0f;
        }
      }
    }
    return vals;
  }

  private void updatePaintColorModel(double alpha) {
    IndexColorModel icm = _cm.getColorModel();
    int n = 256;
    byte[] r = new byte[n];
    byte[] g = new byte[n];
    byte[] b = new byte[n];
    byte[] a = new byte[n];
    icm.getReds(r);
    icm.getGreens(g);
    icm.getBlues(b);
    byte ba = (byte)(255.0*alpha);
    for (int i=0; i<n; ++i) a[i] = ba;
    r[0] = (byte)255;
    b[0] = (byte)255;
    g[0] = (byte)255;
    a[0] = (alpha==1.0)?(byte)255:(byte)0;
    icm = new IndexColorModel(8,n,r,g,b,a);
    _cm.setColorModel(icm);
    _ipg.setColorModel2(icm);
  }

  private static float[][][] loadImage(
    int n1, int n2, int n3, String fileName)
  {
    float[][][] image = new float[n3][n2][n1];
    try {
      ArrayInputStream ais = new ArrayInputStream(fileName);
      ais.readFloats(image);
      ais.close();
      return image;
    } catch (IOException ioe) {
      System.err.println(ioe);
    }
    return image;
  }

  private static EigenTensors3 loadTensors(String fileName) {
    EigenTensors3 et = null;
    try {
      ObjectInputStream ois =
        new ObjectInputStream(new FileInputStream(fileName));
      et = (EigenTensors3)ois.readObject();
      ois.close();
      return et;
    } catch (IOException ioe) {
      System.err.println(ioe);
    } catch (ClassNotFoundException cnfe) {
      System.err.println(cnfe);
    }
    return et;
  }



  public static void main(String[] args) {
    SwingUtilities.invokeLater(new Runnable() {
      public void run() {
        setupSynthetic(); // For a quick synthetic example
        //realData();  // Use on your own data
      }
    });
  }

  private static void realData() {
    String dir = "directory/to/data/";
    String xyz = "image.dat";
    String tensors = "tensors.dat";
    int    n1 = 401,   n2 = 357,   n3 = 161;
    double d1 = 0.004, d2 = 0.025, d3 = 0.025;
    double f1 = 0.600, f2 = 0.000, f3 = 0.000;
    Sampling s1 = new Sampling(n1,d1,f1);
    Sampling s2 = new Sampling(n2,d2,f2);
    Sampling s3 = new Sampling(n3,d3,f3);
    _et = loadTensors(dir+tensors);
    float[][][] image = loadImage(n1,n2,n3,dir+xyz);
    PaintBrush pb = new PaintBrush(s1,s2,s3,_et);
    TensorGuidedPainting tgp = new TensorGuidedPainting(s1,s2,s3,pb,image);
  }
  
  private static void setupSynthetic() {
    int n1 = 101, n2 = 101, n3 = 101;
    float[][][] image =
      FakeData.seismic3d2010A(n1,n2,n3,20.0,10.0,30.0,0.5,0.2);

    // Compute structure tensors
    print("computing structure tensors");
    LocalOrientFilter lof = new LocalOrientFilter(8.0f);
    _et = lof.applyForTensors(image);

    // Compute semblance values
    LocalSemblanceFilter lsf1 = new LocalSemblanceFilter(2,2);
    LocalSemblanceFilter lsf2 = new LocalSemblanceFilter(4,4);
    LocalSemblanceFilter lsf3 = new LocalSemblanceFilter(16,0);

    print("computing 1D semblance");
    float[][][] sm1 =
      lsf1.semblance(LocalSemblanceFilter.Direction3.W,_et,image);
    print("computing 2D semblance");
    float[][][] sm2 =
      lsf2.semblance(LocalSemblanceFilter.Direction3.VW,_et,image);
    print("computing 3D semblance");
    float[][][] sm3 =
      lsf3.semblance(LocalSemblanceFilter.Direction3.UVW,_et,image);

    // Construct metric tensors
    print("computing metric tensors");
      _et.setEigenvalues(sm3,sm2,sm1);

    Sampling s1 = new Sampling(n1);
    Sampling s2 = new Sampling(n2);
    Sampling s3 = new Sampling(n3);
    PaintBrush pb = new PaintBrush(s1,s2,s3,_et);
    TensorGuidedPainting tgp = new TensorGuidedPainting(s1,s2,s3,pb,image);
  }

  // Debugging
  private static void print(Object o) {
    System.out.println(o);
  }

}
