# beam-corrections
Scripts to derive the frequency-dependent beam correction for a source observed in the MeerKAT primary beam or an FBFUSE coherent beam

**Requirements:**
katbeam (https://github.com/ska-sa/katbeam)
mosaic (https://github.com/wchenastro/Mosaic.git)
astropy, matplotlib, numpy

**Usage:**
MeerKAT primary beam correction example:
`python IB-correct.py --band L --freq 1284 --boresight 12:00:00.0 -30:00:00.0 --source 12:01:00.0 -30:00:00.0`
Arguments:  --band specifies UHF, L, or S band receiver
            --freq sets frequency in MHz
            --boresight RA & Dec pointing coordinates
            --source RA & Dec source coordinates
            --beamextent (optional, default 4 deg) sets the size of the IB array. If the source coordinates are further away from the IB centre, make this bigger.

MeerKAT FBFUSE coherent beam correction example:
`python CB-correct.py --r 5 --size 200 --beam 12:00:00 -30:00:00 --source 12:00:00 -30:00:05 --datetime 2015.10.03 15:23:10.000001 --freq 400 --subarray 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63`
Arguments:  --r sets the spatial precision (arcseconds / pixel). Lower value = more precise, but slower.
            --size sets the array size in pixels. You'll get an error if the coordinates don't fit in the array, but bigger = slower.
            --beam RA & Dec pointing centre coordinates
            --source RA & Dec source coordinates
            --datetime specifies date (YYYY.MM.DD) & time (HH:MM:S.) of observation
            --freq specifies observation frequency
            --subarray a comma-separated list of the antennas that were used in the observation
