
import java.awt.BasicStroke;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.FileDialog;
import java.awt.Frame;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.RenderingHints;
import java.awt.Shape;
import java.awt.Stroke;
import java.awt.Window;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;
import java.awt.event.HierarchyEvent;
import java.awt.event.HierarchyListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.geom.AffineTransform;
import java.awt.geom.Line2D;
import java.awt.geom.Path2D;
import java.awt.geom.PathIterator;
import java.awt.geom.Point2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.logging.Logger;

import javax.swing.AbstractButton;
import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.JScrollPane;
import javax.swing.JSpinner;
import javax.swing.JToolBar;
import javax.swing.SpinnerNumberModel;
import javax.swing.SwingUtilities;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.plaf.basic.BasicToolBarUI;


public class Viewer extends JPanel 
    implements MouseMotionListener, MouseListener, 
               ComponentListener, ActionListener {

    private static final Logger logger = 
	Logger.getLogger(Viewer.class.getName());

    static final int THICKNESS = 0;
    static final long LARGE_IMAGE = 500000l;

    static final int SEGMENTS = 1<<0;
    static final int POLYGONS = 1<<1;
    static final int THINNING = 1<<2;
    static final int BITMAP = 1<<3;
    static final int COMPOSITE = 1<<4;
    static final int HISTOGRAM = 1<<5;    
    static final int ALL = SEGMENTS|POLYGONS|THINNING|BITMAP|COMPOSITE|HISTOGRAM;

    static final Color HL_COLOR = new Color (0xdd, 0xdd, 0xdd, 120);
    static Color[] colors = new Color[]{Color.red, Color.blue, Color.black};

    HistogramChart lineHistogram;
    
    FileDialog fileDialog;

    Bitmap bitmap; // original bitmap
    Bitmap thin; // thinned bitmap
    BufferedImage image; // buffered image
    BufferedImage imgbuf; // rendered image

    Collection<Shape> polygons;
    Collection<Path2D> segments;
    
    List<Double> lineLengths;
    List<Line2D> lines;   
    
    Collection<Shape> composites;
    java.util.List<Shape> highlights = new ArrayList<Shape>();

    double cutoffLength=-1;
    
    double sx = 1., sy = 1.;
    AffineTransform afx = new AffineTransform ();

    JPopupMenu popupMenu;

    int show = SEGMENTS|THINNING|BITMAP;
    int available;

    public Viewer ()  {
        addMouseMotionListener (this);
        addMouseListener (this);

        popupMenu = new JPopupMenu ();
        JMenuItem item;
        popupMenu.add(item = new JMenuItem ("Save polygon bitmap"));
        item.setToolTipText("Save highlighted polygon bitmap");
        item.addActionListener(this);
    }

    public Viewer (File file) throws IOException {
        this ();
        load (file);
    }

    public Viewer (File file, double scale) throws IOException {
        this ();
	//sx = (double)(bitmap.width()+3*THICKNESS)/bitmap.width();
	//sy = (double)(bitmap.height()+3*THICKNESS)/bitmap.height();
        load (file, scale);
    }

    public void mouseDragged (MouseEvent e) {
    }
    
    public void mouseMoved (MouseEvent e) {
        if (bitmap == null) {
            return;
        }

        Point pt = e.getPoint();
        highlights.clear();

        if ((show & POLYGONS) == 0) {
            // polygon off
            for (Shape s : polygons) {
                if (Path2D.contains(s.getPathIterator(afx), pt)) {
                    highlights.add(s);
                }
            }
        }

        if (composites != null && (show & COMPOSITE) == 0) {
            // composite off
            for (Shape s : composites) {
                if (Path2D.contains(s.getPathIterator(afx), pt)) {
                    highlights.add(s);
                }
            }
        }

        //logger.info("## "+highlights.size()+" highlights");

        repaint ();
    }

    public void mouseClicked (MouseEvent e) {
        popupGesture (e);
    }
    public void mouseEntered (MouseEvent e) {}
    public void mouseExited (MouseEvent e) {}
    public void mousePressed (MouseEvent e) {
        popupGesture (e);
    }
    public void mouseReleased (MouseEvent e) {
        popupGesture (e);
    }

    boolean popupGesture (MouseEvent e) {
        boolean popup = e.isPopupTrigger();
        if (popup) {
            popupMenu.show(this, e.getX(), e.getY());
        }
        return popup;
    }

    public void actionPerformed (ActionEvent e) {
        String cmd = e.getActionCommand();
        if (cmd.equalsIgnoreCase("save polygon bitmap")) {
            saveHighlightedPolygon ();
        }
    }

    public void componentResized (ComponentEvent e) {
        highlights.clear();
        repaint ();
    }

    public void componentHidden (ComponentEvent e) {}
    public void componentMoved (ComponentEvent e) {}
    public void componentShown (ComponentEvent e) {}

    public void setScale (double scale) {
        if (bitmap != null) {
            sx = scale;
            sy = scale;
            logger.info("scale x: "+sx + " scale y: "+sy);
            setPreferredSize (new Dimension ((int)(sx*bitmap.width()+.5),
                                             (int)(sy*bitmap.height()+.5)));
            resetAndRepaint ();
        }
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

    FileDialog getFileDialog () {
        if (fileDialog == null) {
            fileDialog = new FileDialog 
                ((Frame)SwingUtilities.getAncestorOfClass
                 (Frame.class, this), "Open Image");
        }
        return fileDialog;
    }

    public boolean isAvailable (int flag) {
        return (available & flag) != 0;
    }

    public File load () throws IOException {
        FileDialog fd = getFileDialog ();

        fd.setMode(FileDialog.LOAD);
        fd.setTitle("Open image");
        fd.setVisible(true);
        String name = fd.getFile();
        File file = null;
        if (null != name) {
            load (file = new File (fd.getDirectory(), name));
        }
        return file;
    }

    public void load (File file) throws IOException {
        load (file, Math.min(sx, sy));
    }

    public void load (File file, double scale) throws IOException {
        sx = scale;
        sy = scale;

        available = ALL;

        long start = System.currentTimeMillis();
        try {
            bitmap = Bitmap.readtif(file);
        } 
        catch (Exception e) {
            logger.info("Problem loading file: not valid Tiff? "
                        +"Attempting conversion.");
            bitmap = Bitmap.createBitmap(ImageUtil.decodeAny(file).getData());
        }
        logger.info("## load image "+bitmap.width()+"x"+bitmap.height()+" in "
                    +String.format("%1$.3fs", 
                                   (System.currentTimeMillis()-start)*1e-3));

        start = System.currentTimeMillis();
        polygons = bitmap.polyConnectedComponents();
        logger.info("## generated "+polygons.size()+" connected components in "
                    +String.format("%1$.3fs", 
                                   (System.currentTimeMillis()-start)*1e-3));

        start = System.currentTimeMillis();

        boolean isLarge = polygons.size() > 4000;
        if (!isLarge) {
            composites = new Segmentation (polygons).getComposites();
            logger.info("## generated "+composites.size()+" composites in "
                        +String.format("%1$.3fs", 
                                       (System.currentTimeMillis()-start)*1e-3));
            available |= COMPOSITE;
        }
        else {
            composites.clear();
            available &=~COMPOSITE;
        }

        start = System.currentTimeMillis();
        //thin = bitmap.skeleton();
        thin = bitmap.thin();
        logger.info("## thinning in "
                    +String.format("%1$.3fs", 
                                   (System.currentTimeMillis()-start)*1e-3));

        start = System.currentTimeMillis();
        // segments are generated for thinned bitmap only, since
        //  it can quite noisy on normal bitmap!
        if (!isLarge) {
            segments = thin.segments();
            logger.info("## generated "+segments.size()+" segments in "
                        +String.format("%1$.3fs", 
                                       (System.currentTimeMillis()-start)*1e-3));
            available |= SEGMENTS;
        }
        else {
            available &=~SEGMENTS;
            segments.clear();
        }
        
        lineLengths = new ArrayList<Double>();
        lines= new ArrayList<Line2D>();
        for(Path2D p2:segments){
        	PathIterator pi=p2.getPathIterator(null);
        	double[] prevPt=null;
        	while(!pi.isDone()){
        		
        		double[] coord= new double[2];
        		pi.currentSegment(coord);
        		if(prevPt!=null){
        			Line2D line = new Line2D.Double(coord[0], coord[1], prevPt[0], prevPt[1]);
        			double lineLength=Math.sqrt((coord[0]-prevPt[0])*(coord[0]-prevPt[0])+(coord[1]-prevPt[1])*(coord[1]-prevPt[1]));
        			lineLengths.add(lineLength);
        			lines.add(line);
        			System.out.println(lineLength);
        		}
        		prevPt=coord;
        		pi.next();
        		available |=HISTOGRAM;
        	}
        		
        }
        setPreferredSize (new Dimension ((int)(sx*bitmap.width()+.5),
					 (int)(sy*bitmap.height()+.5)));
        resetAndRepaint ();
        if(lineHistogram==null){
        	final Viewer v = this;
        	this.lineHistogram=new HistogramChart(lineLengths);
        	lineHistogram.addMouseMotionListener(new MouseMotionListener(){
            	@Override
				public void mouseDragged(MouseEvent arg0) {
					// TODO Auto-generated method stub
					
				}
				@Override
				public void mouseMoved(MouseEvent arg0) {
					double cut=v.lineHistogram.getVal(arg0.getX(),arg0.getY());
					System.out.println(v.lineHistogram.getVal(arg0.getX(),arg0.getY()));
					v.cutoffLength=cut;
					v.resetAndRepaint();
				}
            });
        }else{
        	this.lineHistogram.loadData(lineLengths);
        }
        
    }
    

    void resetAndRepaint () {
        imgbuf = null;
        revalidate ();
        repaint ();
    }

    @Override
    protected void paintComponent (Graphics g) {
        if (bitmap == null) {
            return;
        }

	if (imgbuf == null) {
	    imgbuf = ((Graphics2D)g).getDeviceConfiguration()
		.createCompatibleImage(getWidth (), getHeight());
	    Graphics2D g2 = imgbuf.createGraphics();
	    draw (g2);
	    g2.dispose();
	}

        Rectangle r = getBounds ();
        double x = (r.getWidth()-sx*image.getWidth())/2.;
        double y = (r.getHeight()-sy*image.getHeight())/2.;
        Graphics2D g2 = (Graphics2D)g;

	g2.drawImage(imgbuf, (int)(x+.5), (int)(y+.5), null);

        afx.setToTranslation(x, y);
        afx.scale(sx, sy);

        if (!highlights.isEmpty()) {
            g2.setTransform(afx);
            g2.setRenderingHint(RenderingHints.KEY_RENDERING, 
                                RenderingHints.VALUE_RENDER_QUALITY);
            g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, 
                                RenderingHints.VALUE_ANTIALIAS_ON);

            g2.setPaint(HL_COLOR);
            for (Shape s : highlights) {
                g2.fill(s);
            }
        }
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

        if (composites != null && (show & COMPOSITE) != 0) {
            drawComposites (g2);
        }

        if (polygons != null && (show & POLYGONS) != 0) {
            drawPolygons (g2);
        }
        if (this.cutoffLength>0 && (show & HISTOGRAM)!=0){
        	drawColoredLines(g2);
        }
        
        if (segments != null && (show & SEGMENTS) != 0) {
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
        g2.setPaint(HL_COLOR);
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
    void drawColoredLines(Graphics2D g2){
    	Stroke s = g2.getStroke();
		Color rightColor = new Color(HistogramChart.rightColor.getRed(),
				HistogramChart.rightColor.getGreen(),
				HistogramChart.rightColor.getBlue(),
				(int)(HistogramChart.rightColor.getAlpha() * .5));
		Color leftColor = new Color(HistogramChart.leftColor.getRed(),
				HistogramChart.leftColor.getGreen(),
				HistogramChart.leftColor.getBlue(),
				(int)(HistogramChart.leftColor.getAlpha() * .5));
    	g2.setStroke(new BasicStroke(5.0f));
    	for (int j=0;j<lines.size();j++) {
    	    if(lineLengths.get(j)>=cutoffLength){
    	    	g2.setColor(rightColor);
    	    	
    	    }else{
    	    	g2.setColor(leftColor);
    	    }
    	    g2.draw(lines.get(j));
    	}
    	g2.setStroke(s);
    	
    }

    void drawSegments (Graphics2D g2) {
	int i = 0;
	float[] seg = new float[6];
	for (Path2D p : segments) {
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

    void saveHighlightedPolygon () {
        if (highlights.isEmpty()) {
            return;
        }
        
        FileDialog fd = getFileDialog ();
        fd.setMode(FileDialog.SAVE);
        fd.setTitle("Save bitmap polygon as...");
        fd.setVisible(true);
        String name = fd.getFile();
        if (null != name) {
            File file = new File (fd.getDirectory(), name);
            Bitmap poly = bitmap.crop(highlights.iterator().next());
            try {
                poly.write(file);
            }
            catch (IOException ex) {
                JOptionPane.showMessageDialog
                    (this, "Can't save polygon to file \""+file+"\"!",
                     "Error", JOptionPane.ERROR_MESSAGE);
            }
        }
    }
    
    //simple component to display a histogram of values
	static class HistogramChart extends JPanel implements ComponentListener{
		static Color leftColor=Color.GREEN;
		static Color rightColor=Color.MAGENTA;
		static Color defColor = Color.DARK_GRAY;
		
		Collection<Double> _values;
		
		double _max;
		double _min;
		int[] _histogram;
		int _largestFreq;
		int _buckets=50;
		boolean isLog=true;
		BufferedImage imgbuf;
		double _cutoff=-1;
		public HistogramChart() {
			this(null);
		}
		
		public HistogramChart(Collection<Double> values) {
			if(values!=null)
				loadData(values);
			this.addComponentListener(this);
		}
		public void loadData(Collection<Double> values){
			imgbuf=null;
			_values = values;
			processData();
			setPreferredSize (new Dimension ((int)(100),
					 (int)(100)));
		}
		private void processData(){
			_max = Double.MIN_VALUE;
	    	_min = Double.MAX_VALUE;
	    	for(double d:_values){
	    		if(d>_max){
	    			_max=d;
	    		}
	    		if(d<_min){
	    			_min=d;
	    		}
	    	}
	    	double range = _max-_min;
	    	_histogram = new int[_buckets+1];
	    	_largestFreq=0;
	    	for(double d:_values){
	    		int n=++_histogram[(int)(((d-_min)/(_max-_min))*_buckets)];
	    		if(n>_largestFreq){
	    			_largestFreq=n;
	    		}
	    	}
		}
		

		@Override
		protected void paintComponent(Graphics g) {
			
			if (imgbuf == null) {
			    imgbuf = ((Graphics2D)g).getDeviceConfiguration()
				.createCompatibleImage(getWidth (), getHeight());
			    Graphics2D g2 = imgbuf.createGraphics();
			    drawDistribution(g2);
			    g2.dispose();
			}
			Rectangle r = getBounds ();
	        Graphics2D g2 = (Graphics2D)g;
	        g2.drawImage(imgbuf, (int)(0), (int)(0), null);
		}
		void drawDistribution(Graphics2D g2){
			g2.setColor(Color.white);
			g2.fillRect(0, 0, getWidth(), getHeight());
			g2.setColor(defColor);
			
			double width=this.getWidth()/_buckets;
			double bottom=getHeight();
			double maxHeight=getHeight();
			double logfac=maxHeight/Math.log(_largestFreq+1);
			int x=this.getWidth();
			if(_cutoff>0){
				g2.setColor(leftColor);
				x=(int)((_cutoff-_min)*this.getWidth()/(_max-_min));
			}
	    	for(int i=0;i<_histogram.length;i++){
	    		if(i>(int)(((_cutoff-_min)/(_max-_min))*_buckets)-1){
	    			g2.setColor(rightColor);
	    		}
	    		double unitheight;
	    		if(!isLog){
	    			unitheight=(maxHeight*_histogram[i])/_largestFreq;
	    		}else{
	    			unitheight=Math.log(_histogram[i]+1)*logfac;
	    		}
	    		
	    		g2.fillRect((int)(i*width), (int)(bottom-unitheight), (int)width, (int)unitheight+1);
	    	}
	    	
	    	g2.setColor(Color.black);
	    	g2.drawLine(x,0,x,(int)bottom);
	    	g2.drawString((int)(_cutoff)+"",x, (int)(bottom/2));
	    	
	    }
		public double getVal(double x, double y){
			_cutoff=(((_max-_min)*x)/this.getWidth()+_min);
			imgbuf=null;
			repaint();
			return _cutoff;
		}

		@Override
		public void componentHidden(ComponentEvent arg0) {
			// TODO Auto-generated method stub
			
		}

		@Override
		public void componentMoved(ComponentEvent arg0) {
			// TODO Auto-generated method stub
			
		}

		@Override
		public void componentResized(ComponentEvent arg0) {
			// TODO Auto-generated method stub
			imgbuf=null;
			this.revalidate();
			this.repaint();
		}

		@Override
		public void componentShown(ComponentEvent arg0) {
			// TODO Auto-generated method stub
			
		}
	}

    static class ViewerFrame extends JFrame 
        implements ActionListener, ChangeListener {
        Viewer viewer;
        JToolBar toolbar;
        JToolBar toolbar2;

        ViewerFrame (File file, double scale) throws IOException {
            this ();
            setTitle (file.getName());
            viewer.load(file, scale);
           
        }

        ViewerFrame ()  {
            toolbar = new JToolBar ();
            toolbar2= new JToolBar();
            HierarchyListener hl = new HierarchyListener() {
                
                @Override
                public void hierarchyChanged(HierarchyEvent e) {
                    if ((e.getChangeFlags() & HierarchyEvent.PARENT_CHANGED) == 0) return;
                    JToolBar bar = (JToolBar) e.getComponent();
                    if (!((BasicToolBarUI) bar.getUI()).isFloating()) return;
                    final Window topLevel = SwingUtilities.windowForComponent(bar);
                    if (topLevel instanceof JDialog) {
                        ((JDialog) topLevel).setResizable(true);
                    }    
                }
            };
            toolbar2.addHierarchyListener(hl);
            AbstractButton ab;
            toolbar.add(ab = new JButton ("Load"));
            ab.setToolTipText("Load new file");
            ab.addActionListener(this);
            toolbar.addSeparator();

            toolbar.add(ab = new JCheckBox ("Bitmap"));
            ab.putClientProperty("MASK", BITMAP);
            ab.setToolTipText("Show bitmap image");
            ab.setSelected(true);
            ab.addActionListener(this);
            toolbar.addSeparator();

            toolbar.add(ab = new JCheckBox ("Segments"));
            ab.putClientProperty("MASK", SEGMENTS);
            ab.setSelected(true);
            ab.setToolTipText("Show line segments");
            ab.addActionListener(this);

            toolbar.add(ab = new JCheckBox ("Thinning"));
            ab.putClientProperty("MASK", THINNING);
            ab.setSelected(true);
            ab.setToolTipText("Show thinning image");
            ab.addActionListener(this);
            
            toolbar.add(ab = new JCheckBox ("Polygons"));
            ab.putClientProperty("MASK", POLYGONS);
            ab.setToolTipText("Show connected components");
            ab.addActionListener(this);

            toolbar.add(ab = new JCheckBox ("Composites"));
            ab.putClientProperty("MASK", COMPOSITE);
            ab.setToolTipText
                ("Show connected component composites");
            ab.addActionListener(this);
            
            toolbar.add(ab = new JCheckBox ("Histogram"));
            ab.putClientProperty("MASK", HISTOGRAM);
            ab.setToolTipText
                ("Show histogram coloring of lines");
            ab.addActionListener(this);

            toolbar.addSeparator();
            Box hbox = Box.createHorizontalBox();
            hbox.add(new JLabel ("Scale"));
            hbox.add(Box.createHorizontalStrut(5));
            JSpinner spinner = new JSpinner 
                (new SpinnerNumberModel (1., .1, 50., .2));
            spinner.addChangeListener(this);
            hbox.add(spinner);
            hbox.add(Box.createHorizontalGlue());
            toolbar.add(hbox);
            
            
            JPanel pane = new JPanel (new BorderLayout (0, 2));
            pane.add(toolbar, BorderLayout.NORTH);
            pane.add(toolbar2,BorderLayout.SOUTH);
            pane.add(new JScrollPane (viewer = new Viewer ()));
            getContentPane().add(pane);
            pack ();
            setDefaultCloseOperation (JFrame.EXIT_ON_CLOSE);
        }

        public void actionPerformed (ActionEvent e) {
            String cmd = e.getActionCommand();
            AbstractButton ab = (AbstractButton)e.getSource();
            boolean show = ab.isSelected();

            if (cmd.equalsIgnoreCase("load")) {
                File file = null;
                try {
                    file = viewer.load();
                    if (file != null) {
                        setTitle (file.getName());
                        for (Component c : toolbar.getComponents()) {
                            if (c instanceof AbstractButton) {
                                ab = (AbstractButton)c;
                                Integer mask = 
                                    (Integer)ab.getClientProperty("MASK");
                                if (mask != null) {
                                    ab.setEnabled(viewer.isAvailable(mask));
                                }
                            }
                        }
                        repaint ();
                    }
                }
                catch (Exception ex) {
                    ex.printStackTrace();
                    JOptionPane.showMessageDialog
                        (this, "Can't load file \""+file+"\"; perhaps "
                         +"it's not a 1 bpp TIFF image?", "Error", 
                         JOptionPane.ERROR_MESSAGE);
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
            else if (cmd.equalsIgnoreCase("histogram")) {
                viewer.setVisible(HISTOGRAM, show);
            }
        }

        public void stateChanged (ChangeEvent e) {
            JSpinner spinner = (JSpinner)e.getSource();
            viewer.setScale(((Number)spinner.getValue()).doubleValue());
            repaint ();
        }

        public void load (File file, double scale) throws IOException {
            viewer.load(file, scale);
            toolbar2.add(new JLabel("Line Lengths Histogram:"), BorderLayout.NORTH);
            toolbar2.add(viewer.lineHistogram, BorderLayout.NORTH);
            
            repaint ();
        }
    }
    

    static JFrame createApp (String name, double scale) throws IOException {
        logger.info("Loading "+name+"; scale="+scale+"...");
        ViewerFrame vf = new ViewerFrame (new File (name), scale);
	return vf;
    }

    public static void main (final String[] argv) {
	SwingUtilities.invokeLater(new Runnable () {
		public void run () {
		    try {
                        final ViewerFrame vf = new ViewerFrame ();
                        if (argv.length > 0) {
                            try {
                                double scale = 1.;
                                if (argv.length > 1) {
                                    scale = Double.parseDouble(argv[1]);
                                    scale = Math.max(scale, 1.);
                                }
                                vf.load(new File (argv[0]), scale);
                            }
                            catch (NumberFormatException ex) {
                                logger.warning("Bogus scale value: "+argv[1]);
                            }
                        }
			vf.setVisible(true);
		    }
		    catch (Exception ex) {
			ex.printStackTrace();
		    }
		}
	    });
    }
}
