
import java.util.*;
import java.io.*;
import java.awt.*;
import java.awt.image.*;
import java.awt.event.*;
import javax.swing.*;
import java.awt.geom.*;

import java.util.logging.Logger;
import java.util.logging.Level;


public class Viewer extends JPanel {
    private static final Logger logger = 
	Logger.getLogger(Viewer.class.getName());

    static final int THICKNESS = 1;

    BufferedImage image;
    java.util.List<Shape> bboxes;
    java.util.List<GeneralPath> segments;
    double sx, sy;
    BufferedImage imgbuf;

    public Viewer (File file) throws IOException {
	Bitmap bitmap = Bitmap.readtif(file);
	bboxes = bitmap.connectedComponents(Bitmap.Bbox.Polygon);
	segments = bitmap.segments();
	logger.info(bboxes.size() + " connected components; "
		    +segments.size()+" segments!");

	//sx = (double)(bitmap.width()+3*THICKNESS)/bitmap.width();
	//sy = (double)(bitmap.height()+3*THICKNESS)/bitmap.height();
	sx = 1;
	sy = 1;
	logger.info("scale x: "+sx + " scale y: "+sy);
	image = bitmap.createBufferedImage();
	setPreferredSize (new Dimension ((int)(sx*bitmap.width()+.5),
					 (int)(sy*bitmap.height()+.5)));
    }

    Color[] colors = new Color[]{Color.red, Color.blue, Color.black};

    @Override
    protected void paintComponent (Graphics g) {
	if (imgbuf == null) {
	    imgbuf = ((Graphics2D)g).getDeviceConfiguration()
		.createCompatibleImage(getWidth (), getHeight());
	    Graphics2D g2 = imgbuf.createGraphics();
	    draw (g2);
	    g2.dispose();
	}

	g.drawImage(imgbuf, 0, 0, null);
    }

    void draw (Graphics2D g2) {
	g2.setRenderingHint(RenderingHints.KEY_RENDERING, 
			    RenderingHints.VALUE_RENDER_QUALITY);
	g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, 
			    RenderingHints.VALUE_ANTIALIAS_ON);
	g2.setPaint(Color.white);
	g2.fillRect(0, 0, getWidth(), getHeight());
	    
	g2.drawImage(image, THICKNESS, THICKNESS, null);
	//drawBBoxes (g2);
	drawSegments (g2);
    }

    void drawBBoxes (Graphics2D g2) {
	g2.setPaint(Color.red);
	g2.scale(sx, sy);
	for (Shape b : bboxes) {
	    g2.draw(b);
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
    

    static JFrame createApp (String name) throws IOException {
	JFrame f = new JFrame ();
	File file = new File (name);
	f.setTitle(file.getName());
	Viewer v = new Viewer (file);
	JPanel p = new JPanel (new BorderLayout ());
	p.add(new JScrollPane (v));
	f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	f.getContentPane().add(p);
	f.pack();
	return f;
    }

    public static void main (final String[] argv) {
	if (argv.length == 0) {
	    System.err.println("Usage: Viewer FILE");
	    System.exit(1);
	}

	SwingUtilities.invokeLater(new Runnable () {
		public void run () {
		    try {
			JFrame f = createApp (argv[0]);
			f.setVisible(true);
		    }
		    catch (Exception ex) {
			ex.printStackTrace();
		    }
		}
	    });
    }
}
