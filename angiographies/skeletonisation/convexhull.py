"""
Create convex hull from binary numpy image

Running environment requirements: 

    numpy
    scipy
    skimage
    sitk

"""

import numpy as np
import argparse
import time
import scipy
import skimage
from angiographies.utils.iositk import writeSITK, readNIFTIasSITK
from angiographies.utils.formatconversionsitk import numpyToSITK, SITKToNumpy



def flood_fill_hull(image): 
    '''image is a binary np matrix
    returns a binary image with 1 on the voxels belonging to the convex hull, and the convex hull'''   
    start_time = time.time()
    points = np.transpose(np.where(image))
    hull = scipy.spatial.ConvexHull(points)
    deln = scipy.spatial.Delaunay(points[hull.vertices]) 
    idx = np.stack(np.indices(image.shape), axis = -1)
    out_idx = np.nonzero(deln.find_simplex(idx) + 1)
    out_img = np.zeros(image.shape, np.int32)
    out_img[out_idx] = 1
    print("--- %s seconds ---" % (time.time() - start_time))
    return out_img, hull

def convex_hull(image):
    '''image is a binary np matrix
    returns a binary image with true on the voxels belonging to the convex hull'''
    start_time = time.time()
    hull = skimage.morphology.convex_hull_image(image)
    print("--- %s seconds ---" % (time.time() - start_time))
    return hull.astype(int)


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-ifile", help="path to binary image with points", default="", required=True)
    parser.add_argument("-ofile", help="path to output mask", default="", required=True)

    args = parser.parse_args()
    inputf = args.ifile
    outputf = args.ofile

    
    img = readNIFTIasSITK(inputf)
    print(img.GetOrigin())
    print(img.GetSize())
    points = SITKToNumpy(img)
    nidus, _ = flood_fill_hull(points)

    masksitk = numpyToSITK(nidus)
    masksitk.SetOrigin(img.GetOrigin())
    masksitk.SetSpacing(img.GetSpacing())
    print(masksitk.GetOrigin())
    print(masksitk.GetSize())
    writeSITK(masksitk, outputf)



if __name__ == "__main__":
    main()
