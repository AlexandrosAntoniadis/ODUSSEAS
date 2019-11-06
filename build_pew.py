#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  6 12:59:13 2017

@author: aantoniadis
"""

import argparse
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from PyAstronomy import pyasl

plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['axes.spines.right'] = False
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.linewidth'] = 2
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['ytick.major.width'] = 2


def _parser():
    parser = argparse.ArgumentParser(description='Measure the pseudo-EW')
    parser.add_argument('fname', help='Fits file to read (ARES format)')
    parser.add_argument('w1', help='Lower wavelength bound', type=float)
    parser.add_argument('w2', help='Upper wavelength bound', type=float)
    parser.add_argument('-dw', help='minimum wavelength area', type=float, default=0.4)
    parser.add_argument('-p', '--plot', help='Plot fit', action='store_true', default=True)
    return parser.parse_args()


def read_data(fname):
    
    if np.shape(fname)==(2,):
        fname=fname[0]
    flux = fits.getdata(fname)
    hdr = fits.getheader(fname)
    w0, dw, N = hdr['CRVAL1'], hdr['CDELT1'], hdr['NAXIS1']
    wavelength = w0 + dw * np.arange(N)
    
    if round(dw,4) != 0.010:
        cdelt1 = 0.010
        f2 = interp1d(wavelength, flux, kind='linear')
        wavelength = np.arange(wavelength[0], wavelength[-1], cdelt1)
        flux = f2(wavelength)
    
    return wavelength, flux


def cut_data(w, f, w1=None, w2=None):
    if w1 is None:
        w1 = w[0]
    if w2 is None:
        w2 = w[-1]
    idx = (w1 <= w) & (w <= w2)
    return w[idx], f[idx]


def area_between(f, g, dx):
   
    h = abs(g-f)/g
    A = np.trapz(h, dx=dx)
    
    return A


def pseudo_EW(fname, w1, w2, dw=0.4, plot=False):
    wavelength, flux = read_data(fname)
    wavelength, flux = cut_data(wavelength, flux, w1, w2)

    # Find central wavelength
    idx_fmin = np.argmin(flux)
    critical_wl = wavelength[idx_fmin]

    # Work on left side
    wavelength_left, flux_left = cut_data(wavelength, flux, critical_wl-dw, critical_wl)
    idx_left = np.argmax(flux_left)
    wl_max_flux_left = wavelength_left[idx_left]
    fl_max_flux_left = flux_left[idx_left]

    # Work on right side
    wavelength_right, flux_right = cut_data(wavelength, flux, critical_wl, critical_wl+dw)
    idx_right = np.argmax(flux_right)
    wl_max_flux_right = wavelength_right[idx_right]
    fl_max_flux_right = flux_right[idx_right]
    
    # set the area
    x1, x2 = wl_max_flux_left, wl_max_flux_right
    y1, y2 = fl_max_flux_left, fl_max_flux_right
    g = np.polyfit([x1, x2], [y1, y2], deg=1)
    N = len(flux_left[idx_left::]) + len(flux_right[0:idx_right])
    x = np.linspace(x1, x2, N)
    idx = (x1<=wavelength) & (wavelength<=x2)
    f = flux[idx]
    g = np.poly1d(g)(x)
    area = area_between(f, g, dx=x[1]-x[0])*1000
    print(r'Area of line: {:.2f}mÅ'.format(area))

    if plot:
        plt.figure(figsize=(8, 4))
        plt.plot(wavelength, flux)
        plt.plot(x, g)
        plt.fill_between(x, flux[idx], g, where=g>=flux[idx], alpha=0.3)
        plt.axvline(critical_wl, linestyle='--', color='C3')
        plt.axvline(wl_max_flux_left, color='C2')
        plt.axvline(wl_max_flux_right, color='C2')
        plt.xlabel(r"Wavelength [$\AA$]")
        plt.ylabel("Flux")
    
        plt.tight_layout()
        plt.show()
    return area

def main():
    args = _parser()
    
    fname = args.fname    
    w1, w2 = args.w1, args.w2
    dw = args.dw
    plot = args.plot
    pseudo_EW(fname=fname, w1=w1, w2=w2, dw=dw, plot=plot)  
    


if __name__ == '__main__':
    main()
