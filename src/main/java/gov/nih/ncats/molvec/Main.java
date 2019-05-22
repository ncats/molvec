package gov.nih.ncats.molvec;

import gov.nih.ncats.molvec.ui.Viewer;
import org.apache.commons.cli.*;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.nio.file.Files;

/**
 * Created by katzelda on 5/20/19.
 */
public class Main {

    public static void main(String[] args){

        Options options = new Options();

        options.addOption("h", "help", false,"print usage text");

        options.addOption("gui", false, "Run Molvec in GUI mode. file and scale option may be set to preload file");
        options.addOption("scale", true, "scale of image to show in viewer (only valid if gui mode AND file are specified)");

        options.addOption("f", "file", true,"path of image file to process. Supported formats include png, jpeg, tiff.  This option is required if not using -gui");

        options.addOption("o","out",true, "path of output processed mol. Only valid when not using gui mode. If not specified output is sent to STDOUT");

        CommandLineParser parser = new DefaultParser();
        try {
            // parse the command line arguments
            CommandLine commandLine = parser.parse( options, args );
            if(commandLine.hasOption("h")){
                showHelp(options, System.out);
                return;
            }

            if(commandLine.hasOption("gui")){
                //file and scale
                if(commandLine.hasOption("f")){

                    String filepath = commandLine.getOptionValue("f");
                    String scale= "1";
                    if(commandLine.hasOption("scale")){
                        scale = commandLine.getOptionValue("scale");
                    }
                    Viewer.main(new String[]{filepath, scale});
                }else {
                    Viewer.main(new String[0]);
                }
            }else if(commandLine.hasOption("f")){


                    String mol = Molvec.ocr(new File(commandLine.getOptionValue("f")));
                    if(commandLine.hasOption("o")){
                        File outputFile = new File(commandLine.getOptionValue("o"));
                        File parent = outputFile.getParentFile();
                        if(parent !=null){
                            Files.createDirectories(parent.toPath());
                        }

                        try(PrintWriter writer = new PrintWriter(new FileWriter(outputFile))){
                            writer.println(mol);
                        }
                    }else{
                        System.out.println(mol);
                    }
            }else{
                //invalid
                throw new ParseException("gui mode or file not specified");
            }
        }catch(ParseException e){
            System.err.println( "Parsing failed.  Reason: " + e.getMessage() );
            showHelp(options, System.err);
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(-1);
        }
    }

    private static void showHelp(Options options, PrintStream ps){
        HelpFormatter formatter = new HelpFormatter();
        try(PrintWriter writer = new PrintWriter(ps)) {
            formatter.printHelp(writer, HelpFormatter.DEFAULT_WIDTH, "molvec [ -f <path> | -gui] ", null, options, HelpFormatter.DEFAULT_LEFT_PAD, HelpFormatter.DEFAULT_DESC_PAD, null);

        }
    }
}
