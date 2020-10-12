package gov.nih.ncats.molvec;

import gov.nih.ncats.molvec.internal.util.CachedSupplier;
import gov.nih.ncats.molvec.internal.util.ConnectionTable;
import gov.nih.ncats.molvec.internal.util.GeomUtil;

import java.awt.geom.Rectangle2D;
import java.util.Optional;

public class MolvecOptions {
    private double averageBondLength = 1D;
    private boolean center = true;
    private boolean includeSgroups = true;

    public MolvecOptions averageBondLength(double averageBondLength){
        if(averageBondLength <=0){
            throw new IllegalArgumentException("avg bond length must be > 0");
        }
        this.averageBondLength = averageBondLength;
        return this;
    }

    public MolvecOptions center(boolean center){
        this.center = center;
        return this;
    }
    public MolvecOptions includeSgroups(boolean includeSgroups){
        this.includeSgroups = includeSgroups;
        return this;
    }

    protected MolvecResult computeResult(ConnectionTable ct){
        String mol= ct.toMol(averageBondLength, center, includeSgroups);

        return new Result(mol, CachedSupplier.of(()-> ct.getNodes()
                                                    .stream()
                                                    .map(n->n.getPoint())
                                                    .collect(GeomUtil.convexHull())
                                                    .getBounds2D()));
    }

    private static class Result implements MolvecResult{
        private final String mol;
        private final CachedSupplier<Rectangle2D> boundsSupplier;


        public Result(String mol, CachedSupplier<Rectangle2D> boundsSupplier) {
            this.mol = mol;
            this.boundsSupplier = boundsSupplier;
        }

        @Override
        public Optional<String> getMol() {
            return Optional.of(mol);
        }

        @Override
        public Optional<Rectangle2D> getOriginalBoundingBox() {
            return Optional.of(boundsSupplier.get());
        }

        @Override
        public boolean hasError() {
            return false;
        }

        @Override
        public Optional<Throwable> getError() {
            return Optional.empty();
        }


    }
}
