# K2 Light Curves Class
Code snippets and data for Columbia seminar on measuring rotation periods from K2 light curves

Dependencies: astropy, astroML, astroML_addons, gatspy, matplotlib, numpy, scipy, supersmoother, and https://github.com/stephtdouglas/k2spin (the last is not pip installable - you'll need to download it and put it into your python path.)

The data/ directory contains detrended light curves from the K2 Systematics Correnction (K2SC; Aigrain et al. 2016) and K2 Self-flat fielded (K2SFF; Vanderburg & Johnson 2014) codes. 

To run a Lomb-Scargle Periodogram analysis on a set of lightcurves, use the lombscargle.py script. Note that it will not automatically overwrite files. 

    $ python lombscargle.py k2sc_all.lst k2sc_results.csv

This will output some ugly figures. You will then probably want to move all the figures produced into their own subfolder before running the code on a different set of light curves.

The autocorrelation function (ACF) code snippets can be found in acf.py. Some assembly required. 
