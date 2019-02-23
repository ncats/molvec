package tripod.molvec.algo;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.function.Consumer;
import java.util.stream.Collectors;

public class ShellCommandRunner {
	List<String> commandParams;
	File startDir = new File("./");
	
	public static boolean isWindows(){
		boolean isWindows = System.getProperty("os.name")
				  .toLowerCase().startsWith("windows");
		return isWindows;
		
	}
	private static class StreamGobbler implements Runnable {
	    private InputStream inputStream;
	    private Consumer<String> consumer;
	 
	    public StreamGobbler(InputStream inputStream, Consumer<String> consumer) {
	        this.inputStream = inputStream;
	        this.consumer = consumer;
	    }
	 
	    @Override
	    public void run() {
	        new BufferedReader(new InputStreamReader(inputStream)).lines()
	          .forEach(consumer);
	    }
	    
	    public void kill() throws Exception{
	    	inputStream.close();
	    }
	    
	}
	
	
	public Monitor run() throws IOException, InterruptedException{
		ProcessBuilder builder = new ProcessBuilder();
		if (isWindows()) {
		    builder.command(commandParams.toArray(new String[0]));
		} else {
		    builder.command(commandParams.toArray(new String[0]));
		}
		builder.directory(startDir);
		Process process = builder.start();
		
		return new Monitor(process);
	}
	
	
	
	
	public static class Builder{
		List<String> commandParams = new ArrayList<String>();
		File startDir=null;
		
		
		public Builder command(String... cmd){
			commandParams=Arrays.stream(cmd)
								.collect(Collectors.toList());
			return this;
		}
		
		public Builder activeDir(File f){
			this.startDir=f;
			return this;
		}
		
		public Builder activeDir(String s){
			this.startDir=new File(s);
			return this;
		}
		
		
		public ShellCommandRunner build(){
			ShellCommandRunner scr = new ShellCommandRunner();
			scr.commandParams=this.commandParams;
			scr.startDir=startDir;
			return scr;
		}
	}
	
	public static class Monitor{
		
		ExecutorService s1=Executors.newSingleThreadExecutor();
		ExecutorService s2=Executors.newSingleThreadExecutor();
		
		
		Process p;
		
		PrintWriter pw =null;
		
		Consumer<String> onIn=null;
		Consumer<String> onErr=(l)->System.err.println(l);
		Consumer<Integer> onExit=(i)->{
			try{
				this.kill();
			}catch(Exception e){
				e.printStackTrace();
			}
		};
		
		StreamGobbler streamGobblerIn;
		StreamGobbler streamGobblerErr;
		
		private Monitor(Process process){
			this.p=process;
			streamGobblerIn = 
					new StreamGobbler(process.getInputStream(), l->_onInput(l));
			s1.submit(streamGobblerIn);
			streamGobblerErr = 
					new StreamGobbler(process.getErrorStream(), l->_onErr(l));
			s2.submit(streamGobblerErr);
			pw = new PrintWriter(this.p.getOutputStream());
			
			Runnable r = ()->{
				while(this.p.isAlive()){
					try {
						Thread.sleep(10);
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
				this.onExit.accept(p.exitValue());
			};
			new Thread(r).start();
			
		}
		
		public Monitor onKilled(Consumer<Integer> onExit){
			this.onExit=onExit;
			if(!p.isAlive()){
				onExit.accept(p.exitValue());
			}
			return this;
		}
		
		public void kill() throws Exception{
			this.getProcess().destroy();
			streamGobblerIn.kill();
			streamGobblerErr.kill();
			this.s1.shutdown();
			this.s2.shutdown();
		}
		
		private void _onInput(String line){
			if(onIn!=null){
				onIn.accept(line);
			}
		}
		private void _onErr(String line){
			if(onErr!=null){
				onErr.accept(line);
			}
		}
		
		public OutputStream getOut(){
			return p.getOutputStream();
		}
		public InputStream getIn(){
			return p.getInputStream();
		}
		public InputStream getErr(){
			return p.getErrorStream();
		}
		public Process getProcess(){
			return p;
		}
		
		public Monitor onInput(Consumer<String> oni){
			this.onIn=oni;
			return this;
		}
		
		public Monitor onError(Consumer<String> one){
			this.onErr=one;
			return this;
		}
		
		public synchronized Monitor writeLine(String line){
			pw.println(line);
			pw.flush();
			return this;
		}
		
		
		
		
	}
	
	public static void main(String[] args) throws IOException, InterruptedException{
		Monitor m=(new Builder()).activeDir("/home/tyler/workspace/cnsmpo")
		               .command("osra", "-f sdf", "/home/tyler/workspace/molvec/src/test/resources/regressionTest/usanWrongSet1/cas-229975-97-7.png")
		               .build()
		               
		               .run();
		m.onInput(l->{
							System.out.println(l);
		            	   if(l.equals("$$$$")){
		            		   
		            		   try {
								m.kill();
							} catch (Exception e) {
								// TODO Auto-generated catch block
								e.printStackTrace();
							}
		            		   
		            	   }
		               });
		
		
		
		
	}
	
}
