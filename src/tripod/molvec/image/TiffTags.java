package tripod.molvec.image;

public interface TiffTags {
    static final int TAG_IMAGELENGTH = 257;
    static final int TAG_IMAGEWIDTH = 256;
    static final int TAG_RESOLUTIONUNIT = 296;
    static final int   RESOLUTIONUNIT_NONE = 1; // no absolute unit
    static final int   RESOLUTIONUNIT_INCH = 2; // inch
    static final int   RESOLUTIONUNIT_CENT = 3; // centimeter
    static final int TAG_XRESOLUTION = 282;
    static final int TAG_YRESOLUTION = 283;
    static final int TAG_BITSPERSAMPLE = 258;
    static final int TAG_SAMPLESPERPIXEL = 277;
    static final int TAG_PHOTOMETRIC = 262;
    static final int   PHOTOMETRIC_WHITEISZERO = 0;
    static final int   PHOTOMETRIC_BLACKISZERO = 1;
    static final int   PHOTOMETRIC_RGB = 2;
    static final int   PHOTOMETRIC_PALETTE = 3;
    static final int   PHOTOMETRIC_MASK = 4;
    static final int TAG_STRIPOFFSETS = 273;
    static final int TAG_ROWSPERSTRIP = 278;
    static final int TAG_STRIPBYTECOUNTS = 279;
    static final int TAG_COMPRESSION = 259;
    static final int TAG_IMAGEDESCRIPTION = 270;
    static final int TAG_MAXSAMPLEVALUE = 281;
    static final int TAG_MINSAMPLEVALUE = 280;
    static final int TAG_PLANARCONFIGURATION = 284;
    static final int   PLANARCONFIGURATION_CHUNKY = 1;
    static final int   PLANARCONFIGURATION_PLANAR = 2;
}
