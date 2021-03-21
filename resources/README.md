#Resources For Xclean

The molvec `xclean` experimantal cleaning uses some resources to help it determine whether a given cleaned
molecule or image is likely to be an improvement. These resources can include whitelisted chemical structures,
common chemical features, etc. This is a work-in-progress, and some of these resources can be quite large.

To use a whitelist set of inchikeys for analysis, there should be a `ikeys.txt` file in this `resources` folder.
The `ikeys.txt` file will have 1 inchikey per row, and will be used in the `-xclean -inchiscorer` parameters, or
with the java code:

```
ModifiedMolvecPipeline.setInChIDefaultScorer();
```

Which will setup the inchikeys as useful for the default scoring.
