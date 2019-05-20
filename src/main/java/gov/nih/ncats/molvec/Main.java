package gov.nih.ncats.molvec;

import org.apache.commons.cli.*;

/**
 * Created by katzelda on 5/20/19.
 */
public class Main {

    public static void main(String[] args){

        Options options = new Options();

        options.addOption("h", "help", false,"print usage text");

        options.addOption("gui", false, "Run Molvec in GUI mode.");

        options.addOption("file", "path of image file to process. Supported formats include png, jpeg, tiff");

        CommandLineParser parser = new DefaultParser();
        try {
            // parse the command line arguments
            CommandLine line = parser.parse( options, args );
        }catch(ParseException e){
            System.err.println( "Parsing failed.  Reason: " + e.getMessage() );
        }
    }
}
