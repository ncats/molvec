package gov.nih.ncats.molvec;

import java.awt.geom.Rectangle2D;
import java.util.Optional;

/**
 * Object to hold the information of a single molvec image to structure result.
 *
 * @since 0.9.8
 */
public interface MolvecResult {
    /**
     * Returns a molfile of the chemical structure.
     * @return an Optional containing the molfile as String or Empty Optional if there was an error.
     */
    Optional<String> getMolfile();

    /**
     * Returns a SDfile ( structure-data file)  of the chemical structure.
     * @return an Optional containing the SDfile as String or Empty Optional if there was an error.
     */
    default Optional<String> getSDfile(){
        Optional<String> molfile = getMolfile();
        if(!molfile.isPresent()){
            return molfile;
        }
        //according to sdfile spec
        //If the SDfile only contains structures, there can be no blank line between the last "M END"
        //and the $$$$ delimiter line.
        String mol = molfile.get();
        StringBuilder builder = new StringBuilder(mol.length() + 6);
        return Optional.of(builder.append(mol).append(System.lineSeparator())
                            .append("$$$$").toString());
    }
    /**
     * Get the bounding box for the found atoms in the original coordinate
     * system of the processed image.  The returned bounding box may not
     * match the atom coordinates in the return mol file
     * if any transformation to recenter or adjust bond lengths was
     * used to generate the mol file.
     * @return an Optional containing a Rectangle2D or an empty optional if there was an error.
     */
    Optional<Rectangle2D> getOriginalBoundingBox();

    /**
     * If there was an error during computing the molvec Result.
     * @return
     */
    boolean hasError();

    /**
     * Return the Throwable error if there is one; or empty optional if there is no error.
     * @return
     *
     * @see #hasError()
     */
    Optional<Throwable> getError();

    /**
     * Factory method to create a MolvecResult that has the given error.
     * @param t
     * @return
     */
    static MolvecResult createFromError(Throwable t){
        return new ErrorResult(t);
    }
}
