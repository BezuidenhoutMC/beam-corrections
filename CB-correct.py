#!/usr/bin/env python
#Tiaan Bezuidenhout, 2021. For inquiries: bezmc93@gmail.com

import argparse
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import datetime, sys, argparse

import katpoint
from astropy import wcs
from astropy.coordinates import SkyCoord
import astropy.units  as u
from mosaic.beamforming import PsfSim

def makeKatPointAntenna(antennaString):
    antennaKat = []

    for antenna in antennaString:
        antkat = katpoint.Antenna(antenna)
        antennaKat.append(antkat)

    return antennaKat

def parseOptions(parser):
    parser.add_argument('--ants', nargs=1, metavar="file", help='antenna coodinates files')
    parser.add_argument('--resolution', nargs=1, metavar="asec", help='resolution in arcsecond')
    parser.add_argument('--size', nargs=1, metavar="num", help='width in pixels')
    parser.add_argument('--beam', nargs=2, metavar=('RA/Azi', 'DEC/Alt'), help='beam centre RA DEC in format H:M:S D:M:S')
    parser.add_argument('--source', nargs=2, metavar=('RA/Azi', 'DEC/Alt'), help='source position RA DEC in format H:M:S D:M:S')
    parser.add_argument('--freq', nargs=1, metavar=('freq'), help='observation frequency in MHz')
    parser.add_argument('--datetime', nargs=2, metavar=('date', 'time'), help='observation time in format: 03.10.2015 15:23:10.000001')
    parser.add_argument('--subarray', nargs='+', metavar="num", help='list of antennas, saperated by comma')
    parser.add_argument("--weight", action="store_true",
            help='apply weights to individual antenna, attach weight after the item in --subarray, e.g., 0:0.5, 1:0.7, 2:0.5 ')

    args = parser.parse_args()

    resolution = args.resolution[0]
    size = int(args.size[0])**2
    freq = float(args.freq[0])*1e6

    print('res',resolution, 'size',size,'freq',freq)

    if args.ants is not None:
        with open(args.ants[0], 'r') as antFile:
            antennaCoords = antFile.readlines()
    else:
        with open("mosaic/antenna.csv", 'r') as antFile:
            antennaCoords = antFile.readlines()
    if args.resolution is not None:
        resolution = float(args.resolution[0])
    if args.size is not None:
        size = int(args.size[0])**2.
    if args.source is not None:
        sourceCoord = args.source
    else:
        parser.error("no source coordinates specified, try --source RA Dec")
    if args.beam is not None:
        SC = SkyCoord(args.beam[0],args.beam[1], frame='icrs', unit=(u.hourangle, u.deg))
        beamCoord = [SC.ra.deg,SC.dec.deg]
    else:
        parser.error("no beam coordinates specified, try --beam RA Dec")
    if args.datetime is not None:
        observeTime=datetime.datetime.strptime(args.datetime[0] + " "
                + args.datetime[1], '%Y.%m.%d %H:%M:%S.%f')
    else:
        parser.error("no time specified, try --datetime date time")
    if args.weight is True:
        weights = []
    else:
        weights = None
    if args.subarray is not None:
        subarray = []
        arrayString = "".join(args.subarray)
        ant_weights = arrayString.split(",")
        for ant_weight in ant_weights:
            ant_weight_pair = ant_weight.split(':')
            if weights is not None:
                subarray.append(int(ant_weight_pair[0]))
                if len(ant_weight_pair) > 1:
                    complex_weight = complex(ant_weight_pair[1])
                    if complex_weight.imag == 0:
                        weights.append(float(ant_weight_pair[1]))
                    else:
                        weights.append(complex_weight)
                else:
                    weights.append(1.0)
            else:
                subarray.append(int(ant_weight_pair[0]))
    else:
        subarray = []

    paras  = {"antennaCoords": antennaCoords,
        "sourceCoord": sourceCoord,
        "beamCoord": beamCoord,
        "observeTime": observeTime,
        "freq": freq,
        "subarray":subarray,
        "size":size,
        "resolution":resolution,
        "weights":weights}

    return paras

def plot(o,beam):
    thisDpi = 96.
    matplotlib.rcParams.update({'font.size': 8})

    array = beam.psf.image

    center = beam.psf.bore_sight.equatorial

    shape = np.array(array).shape
    step = o.resolution/3600.
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

    SC = SkyCoord(o.sourceCoord[0],o.sourceCoord[1], frame='icrs', unit=(u.hourangle, u.deg))
    SC_d = [SC.ra.deg,SC.dec.deg]

    SC_px = wcs_properties.all_world2pix([SC_d],1)
    plt.contour(array,levels=[0.5],colors='white',linestyles='dashed')
    plt.scatter(shape[0]/2.,shape[0]/2.,c='cyan',marker='o',s=50,zorder=999)
    plt.scatter(SC_px[0,0],SC_px[0,1],c='red',marker='x',s=200,zorder=999)

    Correction = array[int(np.round(SC_px[0,1])),int(np.round(SC_px[0,0]))]
    print("Beam intensity at %s, %s = %f" % (SC.ra.to_string(u.hour),SC.dec.to_string(u.deg),Correction))

    plt.show()

def captureNegetiveNumber():
    for i, arg in enumerate(sys.argv):
          if (arg[0] == '-') and arg[1].isdigit(): sys.argv[i] = ' ' + arg

def main():
    captureNegetiveNumber()

    parser = argparse.ArgumentParser()
    options = parseOptions(parser)
    o = argparse.Namespace(**options)

    if o.subarray != []:
        antennaKat = makeKatPointAntenna([o.antennaCoords[ant] for ant in o.subarray])
    else:
        antennaKat = makeKatPointAntenna(o.antennaCoords)

    psf = PsfSim(antennaKat, o.freq)
    
    try:
        beam = psf.get_beam_shape(o.beamCoord, o.observeTime, o.size, o.resolution, o.weights)
    except:
        print("Error: size likely too small to fit PSF at chosen resolution")
        exit()
    #beam.plot_psf('beam.png', overlap = 0.55, shape_overlay=True, interpolation=True)

    plot(o,beam)

if __name__ == "__main__":
    main()