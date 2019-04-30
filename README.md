#MolVec
NCATS (chemical) ocr
engine that can a way to vectorize
chemical images into Chemical objects preserving the 2D layout as much as 
possible. The code is still very raw in terms of utility. Please forward
questions and/or problems to nguyenda@mail.nih.gov.

#How To Build
   
   This project has a dependency on Chemkit, to install those dependencies automatically using the default CDK implementation, run this command:

   $ bash mavenInstall.sh

   Once the dependencies are installed, you can build the whole project as a jar file with:

   $ mvn clean pacakge

   Or install the project into your maven repository using:

   $ mvn install
   
##Example Usage

    File image = ...
    Chemical chemical = Molvec.ocr(image);
    
    System.out.println( chemical.toMol() );
    
    
##Async Support

  New in 0.8 MolVec supports asynchronous calls
  
    CompleteableFuture<Chemical> future = Molvec.ocrAsync( image);
    Chemical chem = future.get(10, TimeUnit.SECONDS);
  

##UI
  Molvec Comes with a Swing Viewer you can use to step
  through each step of the structure recognition process

![Primitives](sample1.png)
