package gov.nih.ncats.molvec;

import java.awt.geom.Rectangle2D;
import java.util.Optional;

public interface MolvecResult {
    Optional<String> getMol();

    Optional<Rectangle2D> getOriginalBoundingBox();

    boolean hasError();

    Optional<Throwable> getError();


    static MolvecResult createFromError(Throwable t){
        return new ErrorResult(t);
    }
}
