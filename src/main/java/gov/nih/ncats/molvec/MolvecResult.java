package gov.nih.ncats.molvec;

import java.awt.geom.Rectangle2D;
import java.util.Optional;

/**
 * Object to hold the information of a single molvec image to structure result.
 */
public interface MolvecResult {
    /**
     * Returns a molfile of the chemical structure.
     * @return an Optional containing the molfile as String or Empty Optional if there was an error.
     */
    Optional<String> getMol();

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
