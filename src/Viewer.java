
import java.util.*;
import java.io.*;
import java.awt.*;
import java.awt.image.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;
import java.awt.geom.*;

import java.util.logging.Logger;
import java.util.logging.Level;


public class Viewer extends JPanel {
    private static final Logger logger = 
	Logger.getLogger(Viewer.class.getName());

    static final int THICKNESS = 0;

    static final int SEGMENTS = 1<<0;
    static final int POLYGONS = 1<<1;
    static final int THINNING = 1<<2;
    static final int BITMAP = 1<<3;
    static final int COMPOSITE = 1<<4;

    static Color[] colors = new Color[]{Color.red, Color.blue, Color.black};

    Bitmap bitmap; // original bitmap
    Bitmap thin; // thinned bitmap
    BufferedImage image; // buffered image
    BufferedImage imgbuf; // rendered image

    Collection<Shape> polygons;
    Collection<GeneralPath> segments;
    Collection<Shape> composites;

    double sx, sy;
    int show = SEGMENTS|THINNING|BITMAP;

    public Viewer ()  {
    }

    public Viewer (File file) throws IOException {
        load (file);
    }

    public Viewer (File file, double scale) throws IOException {
	//sx = (double)(bitmap.width()+3*THICKNESS)/bitmap.width();
	//sy = (double)(bitmap.height()+3*THICKNESS)/bitmap.height();
        load (file, scale);
    }

    public void setScale (double scale) {
	sx = scale;
	sy = scale;
	logger.info("scale x: "+sx + " scale y: "+sy);
	setPreferredSize (new Dimension ((int)(sx*bitmap.width()+.5),
					 (int)(sy*bitmap.height()+.5)));
        resetAndRepaint ();
    }

    public void setVisible (int flag, boolean visible) {
        if (visible) {
            show |= flag;
        }
        else {
            show &= ~flag;
        }
        resetAndRepaint ();
    }

    public void load (File file) throws IOException {
        load (file, Math.min(sx, sy));
    }

    public void load (File file, double scale) throws IOException {
        sx = scale;
        sy = scale;

        bitmap = Bitmap.readtif(file);
        long start = System.currentTimeMillis();
	polygons = bitmap.polyConnectedComponents();
        logger.info("## generated "+polygons.size()+" connected components in "
                    +String.format("%1$.3fs", 
                                   (System.currentTimeMillis()-start)*1e-3));
        start = System.currentTimeMillis();

        composites = new Segmentation (polygons).getComposites();
        logger.info("## generated "+composites.size()+" composites in "
                    +String.format("%1$.3fs", 
                                   (System.currentTimeMillis()-start)*1e-3));

        thin = bitmap.skeleton();
        start = System.currentTimeMillis();
        // segments are generated for thinned bitmap only, since
        //  it can quite noisy on normal bitmap!
	segments = thin.segments();
	logger.info("## generated "+segments.size()+" segments in "
                    +String.format("%1$.3fs", 
                                   (System.currentTimeMillis()-start)*1e-3));

	setPreferredSize (new Dimension ((int)(sx*bitmap.width()+.5),
					 (int)(sy*bitmap.height()+.5)));
        resetAndRepaint ();
    }

    void resetAndRepaint () {
        imgbuf = null;
        revalidate ();
        repaint ();
    }

    @Override
    protected void paintComponent (Graphics g) {
	if (imgbuf == null) {
	    imgbuf = ((Graphics2D)g).getDeviceConfiguration()
		.createCompatibleImage(getWidth (), getHeight());
	    Graphics2D g2 = imgbuf.createGraphics();
	    draw (g2);
	    g2.dispose();
	}

        Rectangle r = getBounds ();
	g.drawImage(imgbuf, (int)((r.getWidth()-sx*image.getWidth())/2.+.5), 
                    (int)((r.getHeight()-sy*image.getHeight())/2.+.5), null);
    }

    void draw (Graphics2D g2) {
	g2.setColor(Color.white);
	g2.fillRect(0, 0, getWidth(), getHeight());

        g2.scale(sx, sy);
	g2.setRenderingHint(RenderingHints.KEY_RENDERING, 
			    RenderingHints.VALUE_RENDER_QUALITY);
	g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, 
			    RenderingHints.VALUE_ANTIALIAS_ON);

        
        image = (show & THINNING) != 0 
            ? thin.createBufferedImage() 
            : bitmap.createBufferedImage();

        if ((show & BITMAP) != 0) {
            g2.drawImage(image, THICKNESS, THICKNESS, null);
        }

        if ((show & COMPOSITE) != 0) {
            drawComposites (g2);
        }

        if ((show & POLYGONS) != 0) {
            drawPolygons (g2);
        }

        if ((show & SEGMENTS) != 0) {
            drawSegments (g2);
        }
    }

    void drawPolygons (Graphics2D g2) {
	g2.setPaint(Color.red);
	for (Shape a : polygons) {
	    g2.draw(a);
        }
    }

    void drawComposites (Graphics2D g2) {
        g2.setPaint(new Color (0xdd, 0xdd, 0xdd, 120));
	for (Shape s : composites) {
	    g2.fill(s);
	}

        g2.setPaint(Color.blue);
	for (Shape s : composites) {
	    g2.draw(s);
	}
    }

    void drawNearestNeighbor (Graphics2D g2, Collection<Shape> polygons) {
        ArrayList<Line2D> lines = new ArrayList<Line2D>();
	for (Shape a : polygons) {
            double min = Double.MAX_VALUE;
            Point2D[] line = null;
            for (Shape b : polygons) {
                if (a != b) {
                    Point2D[] vertex = GeomUtil.nearestNeighborVertices(a, b);
                    double d = GeomUtil.length(vertex[0], vertex[1]);
                    if (line == null || d < min) {
                        line = vertex;
                        min = d;
                    }
                }
            }
            lines.add(new Line2D.Double(line[0], line[1]));
	}

        // nearest neighbor
        double dist = 0.;
        ArrayList<Double> dists = new ArrayList<Double>();
        for (Line2D l : lines) {
            double d = GeomUtil.length(l.getP1(), l.getP2());
            dists.add(d);
            dist += d;
        }
        dist /= lines.size();
        Collections.sort(dists);

        double med = 0.;
        if (lines.size() % 2 == 0) {
            int i = lines.size()/2;
            med = (dists.get(i) + dists.get(i-1))/2.;
        }
        else {
            med = dists.get(lines.size()/2);
        }
        logger.info("### NN distance: ave="+dist+" med="+med);

        g2.setPaint(Color.red);
        for (Line2D l : lines) {
            if (GeomUtil.length(l.getP1(), l.getP2()) < dist) {
                g2.draw(l);
            }
        }

        g2.setPaint(Color.green);
        for (Line2D l : lines) {
            if (GeomUtil.length(l.getP1(), l.getP2()) < med) {
                g2.draw(l);
            }
        }
    }

    void drawSegments (Graphics2D g2) {
	int i = 0;
	float[] seg = new float[6];
	for (GeneralPath p : segments) {
	    g2.setPaint(colors[i%colors.length]);
	    g2.draw(p);
	    PathIterator pi = p.getPathIterator(null);
	    while (!pi.isDone()) {
		int type = pi.currentSegment(seg);
		switch (type) {
		case PathIterator.SEG_LINETO:
		case PathIterator.SEG_MOVETO:
		    g2.drawOval((int)(seg[0]-2), (int)(seg[1]-2), 4, 4);
		    break;
		}
		pi.next();
	    }
	    ++i;
	}
    }

    static class ViewerFrame extends JFrame 
        implements ActionListener, ChangeListener {
        Viewer viewer;
        JToolBar toolbar;
        FileDialog fileDialog = new FileDialog (this, "Open Image");

        ViewerFrame (File file, double scale) throws IOException {
            setTitle (file.getName());
            toolbar = new JToolBar ();
            AbstractButton ab;
            toolbar.add(ab = new JButton ("Load"));
            ab.setToolTipText("Load new file");
            ab.addActionListener(this);
            toolbar.addSeparator();

            toolbar.add(ab = new JCheckBox ("Bitmap"));
            ab.setToolTipText("Show bitmap image");
            ab.setSelected(true);
            ab.addActionListener(this);
            toolbar.addSeparator();

            toolbar.add(ab = new JCheckBox ("Segments"));
            ab.setSelected(true);
            ab.setToolTipText("Show line segments");
            ab.addActionListener(this);

            toolbar.add(ab = new JCheckBox ("Thinning"));
            ab.setSelected(true);
            ab.setToolTipText("Show thinning image");
            ab.addActionListener(this);
            
            toolbar.add(ab = new JCheckBox ("Polygons"));
            ab.setToolTipText("Show connected components");
            ab.addActionListener(this);

            toolbar.add(ab = new JCheckBox ("Composites"));
            ab.setToolTipText
                ("Show connected component composites");
            ab.addActionListener(this);

            toolbar.addSeparator();
            Box hbox = Box.createHorizontalBox();
            hbox.add(new JLabel ("Scale"));
            hbox.add(Box.createHorizontalStrut(5));
            JSpinner spinner = new JSpinner 
                (new SpinnerNumberModel (1., .1, 5., .2));
            spinner.addChangeListener(this);
            hbox.add(spinner);
            hbox.add(Box.createHorizontalGlue());
            toolbar.add(hbox);
            
            JPanel pane = new JPanel (new BorderLayout (0, 2));
            pane.add(toolbar, BorderLayout.NORTH);
            pane.add(new JScrollPane (viewer = new Viewer (file, scale)));
            getContentPane().add(pane);
            pack ();
            setDefaultCloseOperation (JFrame.EXIT_ON_CLOSE);
        }

        public void actionPerformed (ActionEvent e) {
            String cmd = e.getActionCommand();
            AbstractButton ab = (AbstractButton)e.getSource();
            boolean show = ab.isSelected();

            if (cmd.equalsIgnoreCase("load")) {
                fileDialog.setVisible(true);
                String file = fileDialog.getFile();
                if (null != file) {
                    try {
                        viewer.load(new File 
                                    (fileDialog.getDirectory(), file));
                        setTitle (file);
                        repaint ();
                    }
                    catch (Exception ex) {
                        ex.printStackTrace();
                        JOptionPane.showMessageDialog
                            (this, "Can't load file \""+file+"\"; perhaps "
                             +"it's not a 1 bpp TIFF image?", "Error", 
                             JOptionPane.ERROR_MESSAGE);
                    }
                }
            }
            else if (cmd.equalsIgnoreCase("bitmap")) {
                viewer.setVisible(BITMAP, show);
                for (Component c : toolbar.getComponents()) {
                    if (c instanceof AbstractButton) {
                        if ("thinning".equalsIgnoreCase
                            (((AbstractButton)c).getText())) {
                            c.setEnabled(show);
                        }
                    }
                }
            }
            else if (cmd.equalsIgnoreCase("segments")) {
                viewer.setVisible(SEGMENTS, show);
            }
            else if (cmd.equalsIgnoreCase("thinning")) {
                viewer.setVisible(THINNING, show);
            }
            else if (cmd.equalsIgnoreCase("polygons")) {
                viewer.setVisible(POLYGONS, show);
            }
            else if (cmd.equalsIgnoreCase("composites")) {
                viewer.setVisible(COMPOSITE, show);
            }
        }

        public void stateChanged (ChangeEvent e) {
            JSpinner spinner = (JSpinner)e.getSource();
            viewer.setScale(((Number)spinner.getValue()).doubleValue());
            repaint ();
        }
    }
    

    static JFrame createApp (String name, double scale) throws IOException {
        logger.info("Loading "+name+"; scale="+scale+"...");
        ViewerFrame vf = new ViewerFrame (new File (name), scale);
	return vf;
    }

    public static void main (final String[] argv) {
	if (argv.length == 0) {
	    System.err.println("Usage: Viewer FILE [SCALE]");
	    System.exit(1);
	}

	SwingUtilities.invokeLater(new Runnable () {
		public void run () {
		    try {
                        double scale = 1.;
                        if (argv.length > 1) {
                            try {
                                scale = Double.parseDouble(argv[1]);
                                scale = Math.max(scale, 1.);
                            }
                            catch (NumberFormatException ex) {
                                logger.warning("Bogus scale value: "+argv[1]);
                            }
                        }

			JFrame f = createApp (argv[0], scale);
			f.setVisible(true);
		    }
		    catch (Exception ex) {
			ex.printStackTrace();
		    }
		}
	    });
    }
}
