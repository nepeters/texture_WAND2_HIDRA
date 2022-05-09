For both WAND2 and HIDRA, there are two sets of reduced data, which contain datafiles (.xye or .h5)
Autoreduced; No out-of-plane (eta) slicing
Texture-reduced; Sliced out-of-plane per manuscript

For HIDRA, the data and meta-data are included in the .h5 file, which has a structure as follows;
*.h5
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

Also, for the HIDRA data the following run numbers are associated with the following samples.
1600 - Steel #6
1601 - Steel #3
1602 - Steel #9

For WAND2 the .xye files are output from the Mantid reduction (mantidreduce_WAND2.py). These are contained within zip files that need to be unzipped first

Each python file has libraries required for the run, these include
matplotlib
lmfit
tqdm
pyFAI
scipy
numpy
h5py

There is also a custom set of functions contained inside fit_utils.py that are used by the routines for both WAND2 and HIDRA.

The run_seqMAUD.py script runs MAUD in batch text mode to speed up and provide consistent refinements for each dataset. 

It calls MAUD by the maud_batch.bat script (Windows), which requires an environment variable (MAUD_PATH) to be set, which is the folder which contains MAUD. You could also just set this folder manually inside this .bat file.

Running the reduce_fit scripts for both instruments will create an indiviudal MATLAB script which uses MTEX to load the pole figures, generate the ODF and perform the analyses.