
import java.awt.geom.*;
import java.util.*;
import java.util.logging.Logger;
import java.util.logging.Level;

public class MolTrace {
    private static final Logger logger = 
	Logger.getLogger(MolTrace.class.getName());

    private static final boolean DEBUG;
    static {
	boolean debug = false;
	try {
	     debug = Boolean.getBoolean("moltrace.debug");
	}
	catch (Exception ex) {
	}
	DEBUG = debug;
    }
    
    
}
