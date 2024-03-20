#!/usr/bin/env python
#Tiaan Bezuidenhout, 2021. For inquiries: bezmc93@gmail.com

import argparse
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import sys, argparse

from astropy.io import fits
from scipy.interpolate import griddata

from astropy import wcs
from astropy.coordinates import SkyCoord
import astropy.units  as u

from katbeam import JimBeam


def parseOptions(parser):
    parser.add_argument('--band', required=True, nargs=1, metavar=('band'), help='Receiver (UHF/L/S)')
    parser.add_argument('--freq', required=True, nargs=1, metavar=('freq'), help='Obs frequency (MHz)')
    parser.add_argument('--boresight', required=True, nargs=2, metavar=('RA/Azi', 'DEC/Alt'), help='boresight RA DEC in format H:M:S D:M:S')
    parser.add_argument('--source', required=True, nargs=2, metavar=('RA/Azi', 'DEC/Alt'), help='source position RA DEC in format H:M:S D:M:S')
    parser.add_argument('--beamextent', required=False, nargs=1, metavar=('beamextent'), default=[4], help='extent of beam map in degrees')

    args = parser.parse_args()

    paras  = {"band": args.band[0],
        "boresight": args.boresight,
        "source": args.source,
        "freq": args.freq[0],
        "beamextent": args.beamextent[0],
        }

    return paras


def plot(o,beam):
    thisDpi = 96.
    matplotlib.rcParams.update({'font.size': 8})

    array = beam

    b = SkyCoord(o.boresight[0], o.boresight[1], frame='icrs',unit=(u.hourangle, u.deg))
    center = [b.ra.deg,b.dec.deg]

    shape = np.array(array).shape

    step=o.beamextent/beam.shape[0]  #resolution in deg / pixel
    wcs_properties = wcs.WCS(naxis=2)
    wcs_properties.wcs.crpix = [shape[1]/2.-0.5, shape[0]/2.-0.5]
    wcs_properties.wcs.cdelt = [-step, step]
    wcs_properties.wcs.crval = center
    wcs_properties.wcs.ctype = ["RA---TAN", "DEC--TAN"]

    f = plt.figure(figsize=(800./thisDpi, 600./thisDpi), dpi=thisDpi)
    axis = f.add_subplot(111,aspect='equal', projection=wcs_properties)
    ims = axis.imshow(array, cmap=plt.cm.inferno, vmin=0, vmax=1,
             origin='lower')

    axis.set_xlabel('RA ($^\circ$)')
    axis.set_ylabel('Dec ($^\circ$)')

    f.tight_layout()
    axis.set_aspect('auto')
    ra = axis.coords[0]
    dec = axis.coords[1]
    ra.set_major_formatter('hh:mm:ss')
    ra.set_ticklabel(size=8)
    dec.set_ticklabel(size=8, rotation="vertical")
    dec.set_ticks_position('l')
    ra.set_ticks_position('b')
    ra.set_axislabel("RA", size=8)
    dec.set_major_formatter('dd:mm:ss')
    dec.set_axislabel("DEC", size=8)
    plt.subplots_adjust(left=0.10, bottom=0.10, right=0.98, top=0.96,
            wspace=0, hspace=0)

    SC = SkyCoord(o.source[0],o.source[1], frame='icrs', unit=(u.hourangle, u.deg))
    SC_d = [SC.ra.deg,SC.dec.deg]

    SC_px = wcs_properties.all_world2pix([SC_d],1)
    plt.contour(array,levels=[0.5],colors='white',linestyles='dashed')
    plt.scatter(shape[0]/2.,shape[0]/2.,c='cyan',marker='o',s=50,zorder=999)
    plt.scatter(SC_px[0,0],SC_px[0,1],c='red',marker='x',s=200,zorder=999)
    plt.title('%s-band primary beam at %s MHz'%(o.band,o.freq))

    Correction = array[int(np.round(SC_px[0,1])),int(np.round(SC_px[0,0]))]
    print("Primary beam correction at %s, %s = %f" % (SC.ra.to_string(u.hour),SC.dec.to_string(u.deg),Correction))

    plt.show()

def captureNegetiveNumber():
    for i, arg in enumerate(sys.argv):
          if (arg[0] == '-') and arg[1].isdigit(): sys.argv[i] = ' ' + arg

def makebeam(o):
    if o.band=='UHF':
        beam=JimBeam('MKAT-AA-UHF-JIM-2020')
    elif o.band=='L':
        beam=JimBeam('MKAT-AA-L-JIM-2020')
    elif o.band=='S':
        beam=JimBeam('MKAT-AA-S-JIM-2020')

    margin=np.linspace(-o.beamextent/2.,o.beamextent/2.,128)
    x,y=np.meshgrid(margin,margin)

    beampixels=beam.I(x,y,o.freq[0])
    return beampixels


def main():
    captureNegetiveNumber()

    parser = argparse.ArgumentParser()
    options = parseOptions(parser)
    o = argparse.Namespace(**options)
    if o.band not in ('UHF','L','S'):
        print('band must be UHF, L, or S')
        print('exiting...')
        exit()

    model_prim = makebeam(o)   

    plot(o,model_prim)

if __name__ == "__main__":
    main()