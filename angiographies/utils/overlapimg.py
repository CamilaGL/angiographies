"""
Nidus extraction tools:
    - Morphological closing-opening (adapted from CHENOUNE 2019)
    - Skeleton spider location (adapted from BABIN 2018)
        - sphere-based returns sphere mask
        - convex hull-based returns anchor points to get convex hull with other library (not compatible with this environment)

Running environment requirements: 

    numpy
    sitk
    vtk
    json
    scipy


"""

import numpy as np
import argparse
from angiographies.utils.formatconversionsitk import SITKToNumpy, numpyToSITK
from angiographies.utils.imageprocessing import getLargestConnected

from angiographies.utils.iositk import readNIFTIasSITK, writeSITK

def overlapMask(img, mask):
    '''Overlap two numpy arrays'''
    npimg = SITKToNumpy(img).astype(np.int32)
    npmask = SITKToNumpy(mask).astype(np.int32)
    npnidus = npimg * npmask
    nidus = numpyToSITK(npnidus)
    #nidusitk and img to get nidus in image
    #img - (nidusitk and img)
    #nidus = getLargestConnected(nidus)
    nidus.SetOrigin(img.GetOrigin())
    nidus.SetSpacing(img.GetSpacing())
    nidus.SetDirection(img.GetDirection())
    # extracted = numpyToSITK(np.add(npimg,-(SITKToNumpy(nidus))))
    # extracted.SetOrigin(img.GetOrigin())
    # extracted.SetSpacing(img.GetSpacing())
    # extracted.SetDirection(img.GetDirection())
    return nidus#, extracted
    


def main():

    parser = argparse.ArgumentParser()
    #parser.add_argument("-s", help="seed acquisition method", default="5", required=False)

    parser.add_argument("-ifile", help="path to segmentation", default="", required=True)
    parser.add_argument("-imask", help="path to nidus mask", default="", required=True)
    parser.add_argument("-onidus", help="path to output nidus folder", default="", required=True)
    parser.add_argument("-ovessels", help="path to output nifti with the rest", default=None, required=False)

    args = parser.parse_args()
    inputf = args.ifile
    imask = args.imask
    onidus = args.onidus
    ovesselsf = args.ovessels

    img = readNIFTIasSITK(inputf)
    mask = readNIFTIasSITK(imask)
    nidus = overlapMask(img, mask)
    
    #TODO add nidus folder and nidus name
    writeSITK(nidus,onidus)
    # if ovesselsf is not None:
    #     writeSITK(ovesselsf,vessels)



if __name__ == "__main__":
    main()
