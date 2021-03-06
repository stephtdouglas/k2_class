from __future__ import print_function, division

import sys
import logging

import matplotlib
matplotlib.use("agg")
logging.basicConfig(level=logging.WARNING)
from scipy import fftpack
from astropy.convolution import convolve as ap_convolve
from astropy.convolution import Box1DKernel, Gaussian1DKernel
import astropy.io.ascii as at
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import argrelextrema

from kep_io import k2sc_io

def run_fft(times,yvals,modes=np.arange(1,15,3),input_period=None,plot=False):
    """
    runs a FFT and plots the results

    times should be in seconds

    (shouldn't need this for Columbia class)

    """

    plot_ymin,plot_ymax = np.percentile(yvals,[1,99])

    y_fft = fftpack.fft(yvals)
    n = len(y_fft)
    # Output has 0th value is the zero-frequency term,
    # [1:n/2+1] has positive frequency terms where n is the length of the output array
    # [n/2+1:] has negative frequency terms

    # Transform the FFT results into real frequency space (i.e. use the real time spacing),
    day_to_sec = 24*3600.0
    delt = 1.0/(len(times)*np.average(np.diff(times))*day_to_sec)

    freq = np.zeros(n)
    #print n,n/2+1
    if (n%2)==1:
        max_freq = n/2+1
        #print max_freq
        freq[1:max_freq] = np.arange(1,max_freq) * delt
        freq[max_freq:] = np.arange(1,max_freq)[::-1] * delt
    else:
        max_freq = n/2
        #print max_freq
        freq[1:max_freq] = np.arange(1,max_freq) * delt
        freq[max_freq:] = np.arange(1,max_freq+1)[::-1] * delt

    #then transform to period space.
    period = 1.0/freq/day_to_sec # in days
    y_fft1,period1 = y_fft[1:],period[1:]
    # Find its peak, and note it on the plot.
    fund_loc = np.argmax(abs(y_fft1))
    fund_period = period1[fund_loc]
    print_period = "Prot = {:.2f}".format(fund_period)


    if plot:
        plt.figure(figsize=(12,8))

        # Plot the subset of the data,
        # along with fourier series truncated after different numbers of modes

        ax1 = plt.subplot(221)
        ax1.plot(times,yvals,"k.")
        ax1.set_xlabel("time (d)")
        ax1.set_ylabel("flux ")
        for k in np.arange(1,15,3):

            y_fft2 = np.copy(y_fft)
            y_fft2[k + 1:-k] = 0

            y_fit = fftpack.ifft(y_fft2).real

            ax1.plot(times,y_fit,lw=2,label=str(k))
        ax1.legend(loc=3,borderaxespad=0.0,title="# of modes",ncol=5,
             handletextpad=0.05,columnspacing=0.25)
        ax1.set_ylim(plot_ymin*0.98,plot_ymax)


        # Plot the power spectrum as a function of period and note the peak
        ax2 = plt.subplot(222)
        ax2.plot(period,abs(y_fft))
        ax2.set_xlabel("Period (d)")
        ax2.set_ylabel("PSD")
        ax2.plot((fund_period,fund_period),ax2.get_ylim(),'r--',
                 label=print_period)
        ax2.tick_params(labelleft=False,labelright=True)
        if input_period:
            ax2.plot((input_period,input_period),ax2.get_ylim(),"g:",lw=2,
                     label="Input={:.2f}".format(input_period))
        ax2.legend()

        # Phase-fold the light-curve and plot the result
        phase = times/fund_period - np.asarray((times/fund_period),np.int)

        ax3 = plt.subplot(223)
        ax3.plot(phase,yvals,'k.')
        ax3.set_xlabel("Phase")
        ax3.set_ylabel("Flux")
        ax3.set_ylim(plot_ymin*0.98,plot_ymax)
        #print plot_ymin,plot_ymax

    return fund_period

def acf(times,yvals):
    """
    Compute the autocorrelation function for an evenly-sampled time-series

    Inputs:
    -------
    times: array of epochs/time points

    yvals: array of fluxes corresponding to times

    Returns:
    --------
    periods: array of lag times

    ACF: autocorrelation function corresponding to the lags

    """

    # And the number of observations
    # (the maximum lag that will be tested is half the number of observations)
    N = len(yvals)
    max_lag = N//2
    lags = np.arange(max_lag)

    # Calculate values for normalizing the ACF
    median_yval = np.median(yvals)
    norm_term = np.sum((yvals - median_yval)**2)
    #print median_yval,norm_term,max_lag

    # Compute the ACF following McQuillan+2013
    ACF0 = [np.sum((yvals[:N-j] - median_yval)*(yvals[j:] - median_yval)) for j in lags]
    ACF1 = ACF0/norm_term

    # smooth the ACF
    gauss_kernel = Gaussian1DKernel(18,x_size=55)
    ACF = ap_convolve(ACF1, gauss_kernel,boundary="extend")
    #ACF = ACF1

    # Compute the lags in index-space to time-space
    cadence = np.median(np.diff(times))
    periods = cadence*lags

    return periods,ACF


def find_prot(periods,ACF):
    """
    Determines the Prot from an ACF, using procedure in McQuillan et al. (2013)

    Inputs:
    -------
    periods: array of lag times

    ACF: autocorrelation function corresponding to the lags

    Returns:
    --------
    best_period:
         the period corresponding to the first local maximum in the ACF,
         if there are no peaks in the ACF, returns -1

    """

    # Find all local maxima in the ACF. If none, return -1

    max_loc = argrelextrema(ACF,np.greater,order=5)
    #print "max_loc",max_loc
    #print "edge",len(periods)
    if len(max_loc)==0:
        return -1
    max_per = periods[max_loc[0]]
    #print "max_per",max_per
    max_ACF = ACF[max_loc[0]]
    #print "max_acf",max_ACF

    # Find all local minima in the ACF.

    min_loc = argrelextrema(ACF,np.less,order=5)
    #print "min_loc",min_loc
    min_per = periods[min_loc[0]]
    #print "min_per",min_per
    min_ACF = ACF[min_loc[0]]
    #print "min_acf",min_ACF

    ### Find peak heights
    ## Ignore first peak if it's close to 0
    if min_per[0]<1:
        peak_heights = np.zeros(len(max_per)-1)
        per_with_heights = max_per[1:]
        max_ACF_with_heights = max_ACF[1:]
    else:
        peak_heights = np.zeros(len(max_per))
        per_with_heights = max_per
        max_ACF_with_heights = max_ACF

    ## Ignore last peak if there are no minima to the right of it
    while len(np.where(min_per>per_with_heights[-1])[0])==0:
        peak_heights = peak_heights[:-1]
        per_with_heights = per_with_heights[:-1]
        max_ACF_with_heights = max_ACF_with_heights[:-1]
        if len(peak_heights)==0:
            print("No local minima to the right of any local maxima")
            return -1

    for i,max_p in enumerate(per_with_heights):
        # find the local minimum directly to the left of this maximum
        min_left = np.where(min_per<max_p)[0]
        min_loc_1 = min_left[-1]

        # find the local minimum directly to the right of this maximum
        min_right = np.where(min_per>max_p)[0]
        min_loc_2 = min_right[0]
        #print min_per[min_loc_1],max_p,min_per[min_loc_2]
        height1 = max_ACF_with_heights[i] - min_ACF[min_loc_1]
        height2 = max_ACF_with_heights[i] - min_ACF[min_loc_2]
        peak_heights[i] = (height1 + height2)/2.0
    #print peak_heights

    if (len(peak_heights)>1) and (peak_heights[1]>peak_heights[0]):
        # if the second peak is highest, the first peak is probably
        # a half-period alias, so take the second peak.
        best_period = per_with_heights[1]
    else:
        # if the first peak is highest, it's most likely the period
        best_period = per_with_heights[0]

    return best_period


def run_acf(times,yvals,input_period=None,plot=False,output_filename="acf.png"):
    """
    Run the ACF on the input light curve, and plot the result if desired

    Inputs:
    -------
    times: array of epochs/time points

    yvals: array of fluxes corresponding to times

    input_period: (optional) known/existing period for comparison

    plot: (optional) True/False whether to plot outputs

    output_filename: (required if plot=True) filename for resulting plot

    Returns:
    --------
    peak_loc:
         the period corresponding to the first local maximum in the ACF,
         if there are no peaks in the ACF, returns -1


    """

    plot_ymin,plot_ymax = np.percentile(yvals,[1,99])

    periods, ACF = acf(times,yvals)
#    # find the maximum of the first peak
#    peak_locs = argrelextrema(ACF,np.greater,order=5)
#    #print periods[peak_locs[0]]
    peak_loc = find_prot(periods,ACF)
    print_period = "Prot = {:.2f}".format(peak_loc)

    if plot:
        plt.figure(figsize=(10,8))

        ax1 = plt.subplot(221)
        ax1.plot(times,yvals,'k-')
        ax1.set_ylabel("normalized flux",fontsize="large")
        ax1.set_xlabel("Time (d)",fontsize="large")

        ax2 = plt.subplot(222)
        ax2.plot(periods,ACF)
        ax2.set_xlabel(r"$\tau_K$",fontsize="x-large")
        ax2.set_ylabel("ACF",fontsize="large")
        plot2_ymin,plot2_ymax = ax2.get_ylim()
        ax2.set_ylim(plot2_ymin,plot2_ymax)
        if input_period:
            ax2.plot((input_period,input_period),(plot2_ymin,plot2_ymax),"g:",lw=2,label="Input={:.2f}".format(input_period))

        ax2.plot((peak_loc,peak_loc),(plot2_ymin,plot2_ymax),'r--',label=print_period)
        #ax2.plot(periods[peak_locs[0]],ACF[peak_locs[0]],'ro')
        ax2.legend()
        ax2.tick_params(labelleft=False,labelright=True)

        # Phase-fold the light-curve and plot the result
        phase = times/peak_loc - np.asarray((times/peak_loc),np.int)

        ax3 = plt.subplot(223)
        ax3.plot(phase,yvals,'k.')
        ax3.set_xlabel("Phase")
        ax3.set_ylabel("Flux")
        ax3.set_ylim(plot_ymin*0.98,plot_ymax)
        #print plot_ymin,plot_ymax

        plt.savefig(output_filename)

    return peak_loc

if __name__=="__main__":

    listfile = at.read(sys.argv[1])
    file_list = listfile["filename"]

    for filename in file_list:
        print(filename)
        t, f = k2sc_io(filename)
        run_acf(t,f,plot=True,output_filename=filename.split("/")[-1].replace(".fits",".png"))
