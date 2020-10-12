package gov.nih.ncats.molvec;

import java.awt.geom.Rectangle2D;
import java.util.Objects;
import java.util.Optional;

class ErrorResult implements MolvecResult{

    private final Throwable t;

    public ErrorResult(Throwable t) {
        this.t = Objects.requireNonNull(t);
    }

    @Override
    public Optional<String> getMol() {
        return Optional.empty();
    }

    @Override
    public Optional<Rectangle2D> getOriginalBoundingBox() {
        return Optional.empty();
    }

    @Override
    public boolean hasError() {
        return true;
    }

    @Override
    public Optional<Throwable> getError() {
        return Optional.of(t);
    }
}
