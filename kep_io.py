from __future__ import print_function, division
from datetime import date
import logging
import pickle

logging.basicConfig(level=logging.WARNING)

import numpy as np
from astropy.io import fits
import astropy.io.ascii as at


def k2sff_io(filename, ext):
    """ Read in a K2SFF light curve file, and return the time and flux
    arrays from the selected extension.

    Removes timepoints where the MOVING flag is set, indicating data taken
    during a thurster fire.

    Inputs
    ------
    filename: string
        a valid K2SFF file, including full or relative path

    ext: int or string
        The light curve extension to return. If 1, 0, or "best", will return
        the BESTAPER extension, which corresponds to the best light curve as
        determined by Vanderberg et al.

    Returns:
    --------
    time, flux: arrays

    """
    hdu = fits.open(filename)
    #print(hdu.info())
    #print(hdu[ext].data.dtype)

    if ext=="best":
        ext=1
    elif ext==0:
        logging.warning("There is no light curve in extension 0; \n"
                        "returning best aperture")
        ext=1

    table = hdu[ext].data
    time = table["T"][table["MOVING"]==0]
    flux = table["FCOR"][table["MOVING"]==0]
    hdu.close()
    return time,flux

def k2sc_io(filename):
    """ Read in a K2SC light curve file, and return the time and flux.

    Flux is calculated as the completely detrended flux plus the
    time-dependent trend (i.e., only the position-dependent trend is removed),
    and then normalized.

    Inputs
    ------
    filename: string
        a valid K2SC file, including full or relative path

    Returns:
    --------
    time, flux: arrays

    """
    trend_t, trend_p, flux, time = "trtime", "trposi", "flux", "time"

    with fits.open(filename) as hdu:
        #print(hdu.info())

        table = hdu[1].data
        #print(table.dtype)
        good = np.isfinite(table[flux]) & (np.isfinite(table[trend_t]))
        med_trend = np.median(table[trend_t][good])
        time = table[time][good]
        flux = table[flux][good] + table[trend_t][good] - med_trend

    return time,flux
