# texture_WAND2_HIDRA
WAND2 and HIDRA texture analysis codes <br />

### Required python libraries

matplotlib <br />
lmfit <br />
tqdm <br />
pyFAI <br />
scipy <br />
numpy <br />
h5py <br />

There is also a custom set of functions contained inside fit_utils.py that are used by the routines for both WAND2 and HIDRA.

### Data

For both WAND2 and HIDRA, there are two sets of reduced data, which contain datafiles (.xye or .h5) <br />
Auto-reduced: No out-of-plane (eta) slicing <br />
Texture-reduced: Sliced out-of-plane per manuscript <br />

### HIDRA

For HIDRA, the data and meta-data are included in the .h5 file, which generally has a structure as follows; <br />

    *.h5 <br />
    ├── instrument 
    │   ├── calibration 
    │   ├── monochomator setting 
    │   │   └── wave length ** 
    │   ├── efficiency calibration 
    │   └── geometry setup 
    ├── raw data 
    │   ├── sub-runs 
    │   └── logs 
    │       ├── chi    ** 
    │       ├── phi    ** 
    │       └── 2thetaSetpoint ** 
    ├── mask 
    ├── reduced diffraction data 
    │   ├── eta_-5.0 ** 
    │   ├── eta_0.0  ** 
    │   ├── eta_5.0  ** 
    │   └── 2theta   ** 
    └── peaks 
    
Also, for the HIDRA data the following run numbers are associated with the following samples. <br />
1600 - Steel #6 <br />
1601 - Steel #3 <br />
1602 - Steel #9 <br />

### WAND2

For WAND2 the .xye files are output from the Mantid reduction (mantidreduce_WAND2.py). These are contained within zip files that need to be unzipped first

The run_seqMAUD.py script runs MAUD in batch text mode to speed up and provide consistent refinements for each dataset. It calls MAUD by the maud_batch.bat script (Windows), which requires an environment variable (MAUD_PATH) to be set, which is the folder which contains MAUD. You could also just set this folder manually inside this .bat file. Running the reduce_fit scripts for both instruments will create an indiviudal MATLAB script which uses MTEX to load the pole figures, generate the ODF and perform the analyses.

Check with the instrument scientist for the most up to date procedure, these are simply meant as a supplement to a forthcoming publication in JACr.
