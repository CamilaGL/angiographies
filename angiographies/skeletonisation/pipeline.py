"""
Process angiography segmentation to perform nidus extraction.
vmtknetworkextraction only works with vtk < 9.1 or a huge amount of RAM.

Running environment requirements: 

    numpy
    scipy
    skimage
    sitk

"""
import os
import argparse
import time
from angiographies.utils.iositk import writeSITK, readNIFTIasSITK
from angiographies.utils.iovtk import writeVTKPolydataasVTP, readNIFTIasVTK, readVTPPolydata
from angiographies.utils.formatconversionsitk import numpyToSITK, SITKToNumpy
from angiographies.utils.formatconversionvtk import vtkToNumpy, NumpyToVTK
from angiographies.utils.imageprocessing import gaussianSmoothDiscrete
from angiographies.skeletonisation.orderedthinning2 import binarySegmToBinarySkeleton3
from angiographies.skeletonisation.skeletongraph import binSketoSke2
from angiographies.skeletonisation.polydatamerger import skeToPolyline
from angiographies.skeletonisation.nidusextractor import avmMask
from angiographies.skeletonisation.nidusextractor import extractNidusSphere
from angiographies.skeletonisation.vmtknetwork import getvmtkNetwork
from angiographies.skeletonisation.networkediting import toUniquePointID



def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-ipath", help="path to pipeline folder", default="", required=True)
    parser.add_argument("-case", help="case number", default="", required=True)
    parser.add_argument("-segm", help="which segmentation folder are we navigating", default="", required=True)
    parser.add_argument("-method", help="which skeletonisation are we using (skeleton, vmtk, morphological, all)", default="", required=True)
    parser.add_argument("-spider", help="which nidus extraction method based on spiders are we using (boundingbox, spheres, hull, all)", default="", required=True)
    #parser.add_argument("--gaussian", help="perform gaussian smoothing", action="store_true", required=False, default=False)
    parser.add_argument("--overwrite", help="if exists, overwrite skeleton", action="store_true", required=False, default=False)
    parser.add_argument("--grayscale", help="perform grayscale ordered thinning", action="store_true", required=False, default=False)
    parser.add_argument("-rad", help="sphere radius for morphological extraction", type=int, default=15, required=False)
    parser.add_argument("-suffix", help="add suffix to automatic filename", default=None, required=False)
    parser.add_argument("-next", help="discard the n first spiders with max degree", type=int, default=0, required=False)

    args = parser.parse_args()
    ipath = args.ipath
    case = args.case
    segm = args.segm
    method = args.method
    spider = args.spider
    #gaussian = args.gaussian
    grayscale = args.grayscale
    overwrite = args.overwrite
    rad = args.rad
    suffix = args.suffix
    next = args.next

    origin = None
    spacing = None
    shape = None
    methods = ["skeleton", "vmtk", "morphological"] if method == "all" else [method]
    niduses = ["boundingbox", "spheres", "hull"] if spider == "all" else [spider]
    thinmethod = "edt" if not grayscale else "gray"

    segmfile = os.path.join(ipath, "segmentation", segm, case+".nii.gz")

    if os.path.isfile(segmfile): #if we don't have a segmentation file, we're doing nothing

        print("Reading case", segmfile)
        start_time = time.time()

        outputpath = os.path.join(ipath, case) #we're saving everything inside a folder with this case's name
        print("output path", outputpath)
        os.makedirs(outputpath, exist_ok=True) #if the folder doesn't exist, we make it

        filename = case + "_" + segm #this is how our filename starts

        for method in methods: #loop through all methods 

            if method == "vmtk" or method == "skeleton": #we're doing spiders
                polydata = None
                origin = (0,0,0) #origin is the same for both (vmtk loses origin and skeleton doesn't have real world coordinates)

                if method == "vmtk":
                    print ("Using vmtk extraction")
                    filenameske = filename+"_"+method #this is our filename for the skeleton
                    print("Reading file here", os.path.join(outputpath, filenameske+".vtp"))
                    img = readNIFTIasVTK(segmfile) #read segmentation to get volume information. TODO: json
                    shape = img.GetDimensions()[::-1]
                    spacing = img.GetSpacing()
                    if not overwrite and os.path.isfile(os.path.join(outputpath, filenameske+".vtp")): #we don't overwrite if skeleton exists
                        print ("Skeleton exists")
                        polydata = readVTPPolydata(os.path.join(outputpath, filenameske+".vtp"))
                    else:
                        print("Vmtk network should be computed independently")
                        return
                        #network = getvmtkNetwork(img)
                        #polydata = toUniquePointID(network)
                        #writeVTKPolydataasVTP(polydata, os.path.join(outputpath, filenameske+".vtp"))
                    
                else:
                    print ("Doing skeleton extraction", thinmethod)
                    filenameske = filename+"_"+method+thinmethod #this is our filename for the skeleton
                    print("Going to save skeleton here", os.path.join(outputpath, filenameske+".vtp"))
                    img = readNIFTIasSITK(segmfile) #read segmentation to get volume information. TODO: json
                    shape = img.GetSize()[::-1]
                    spacing = (1,1,1)
                    if not overwrite and os.path.isfile(os.path.join(outputpath, filenameske+".vtp")): #we don't overwrite if skeleton exists
                        print ("Skeleton exists")
                        polydata = readVTPPolydata(os.path.join(outputpath, filenameske+".vtp"))
                    else:
                        npimgorig=None
                        if grayscale:
                            if os.path.isfile(os.path.join(ipath, "raw", case+"_0000.nii.gz")): #check if grayscale file exists
                                npimgorig = SITKToNumpy(readNIFTIasSITK(os.path.join(ipath, "raw", case+"_0000.nii.gz"))) if grayscale else None
                            else:
                                print("Case does not have grayscale volume") 
                                return                 
                        
                        img = gaussianSmoothDiscrete(img) #perform gaussian before thinning
                        npimg = SITKToNumpy(img)
                        npthinned = binarySegmToBinarySkeleton3(npimg, npimgorig, 0.25)            
                        ske, _ = binSketoSke2(npthinned, grayscale)
                        polydata = skeToPolyline(ske)
                        writeVTKPolydataasVTP(polydata, os.path.join(outputpath, filenameske+".vtp"))

                for nidus in niduses:
                    if nidus == "spheres" or nidus == "boundingbox" or nidus == "hull":
                        if polydata is not None:
                            print ("Doing spiders with", nidus)
                            
                            filenamemask = filenameske+"_"+nidus
                            if suffix is not None:
                                filenamemask = filenamemask+"_"+suffix
                            if os.path.isfile(os.path.join(outputpath, filenamemask+".nii.gz")) and not overwrite:
                                print("Won't overwrite mask")
                            else:
                                print("Going to save mask here", os.path.join(outputpath, filenamemask+".nii.gz"))
                                mask = avmMask(shape, origin, spacing, polydata, nidus, next) 
                                if mask is not None:
                                    masksitk = numpyToSITK(mask)
                                    masksitk.SetOrigin(img.GetOrigin())
                                    masksitk.SetSpacing(img.GetSpacing())
                                    writeSITK(masksitk, os.path.join(outputpath, filenamemask+".nii.gz"))
                                else:
                                    print("Couldn't find nidus")
                        else:
                            print("No polydata")
                    else:
                        print("Invalid spider method")

            elif method == "morphological": #we don't need to skeletonise for this

                print ("Extracting nidus with morphological operations")
                filenamenidus = filename+"_"+method+"_"+str(rad)
                if suffix is not None:
                    filenamenidus = filenamenidus+"_"+suffix
                if os.path.isfile(os.path.join(outputpath, filenamenidus+".nii.gz")) and not overwrite:
                    print("Won't overwrite mask")
                else:
                    img = readNIFTIasSITK(segmfile)
                    masksitk = extractNidusSphere(img, rad)#spheres
                    masksitk.SetOrigin(img.GetOrigin())
                    masksitk.SetSpacing(img.GetSpacing())
                    writeSITK(masksitk, os.path.join(outputpath, filenamenidus+".nii.gz"))

            else:
                print ("Invalid method")
        print("--- %s seconds ---" % (time.time() - start_time))
    else:
        print("Invalid case")







# def main():

#     parser = argparse.ArgumentParser()
#     parser.add_argument("-ipath", help="path to pipeline folder", default="", required=True)
#     parser.add_argument("-case", help="case number", default="", required=True)
#     parser.add_argument("-segm", help="which segmentation folder are we navigating", default="", required=True)
#     parser.add_argument("-method", help="which skeletonisation are we using (skeleton, vmtk, all, none)", default="", required=True)
#     parser.add_argument("-nidus", help="which nidus extraction method are we using (morphological, boundingbox, spheres, hull, all)", default="", required=True)
#     #parser.add_argument("--gaussian", help="perform gaussian smoothing", action="store_true", required=False, default=False)
#     parser.add_argument("--grayscale", help="perform grayscale ordered thinning", action="store_true", required=False, default=False)
#     parser.add_argument("-rad", help="sphere radius for morphological extraction", type=int, default=30, required=False)

#     args = parser.parse_args()
#     ipath = args.ipath
#     case = args.case
#     segm = args.segm
#     method = args.method
#     nidus = args.nidus
#     #gaussian = args.gaussian
#     grayscale = args.grayscale
#     rad = args.rad

#     origin = None
#     spacing = None
#     shape = None
#     skeletons = []
#     niduses = []
#     skeletons = ["skeleton", "vmtk", "none"] if method == "all" else [method]
#     niduses = ["morphological", "boundingbox", "spheres", "hull"] if nidus == "all" else [nidus]
#     thinmethod = "edt" if not grayscale else "gray"
#     polydata = None

#     segmfile = os.path.join(ipath, "segmentation", segm, case+".nii.gz")
    
#     start_time = time.time()
#     for nidus in niduses:
#         if nidus == "morphological": #we don't need to skeletonise for this

#             print ("Extracting nidus with morphological operations")
#             img = readNIFTIasSITK(segmfile)
#             masksitk = extractNidusSphere(img, rad)#spheres
#             outputpath = os.path.join(ipath, "classification", segm, nidus)

#             masksitk.SetOrigin(img.GetOrigin())
#             masksitk.SetSpacing(img.GetSpacing())
#             os.makedirs(outputpath, exist_ok=True)
#             writeSITK(masksitk, os.path.join(outputpath, case+".nii.gz"))

#         elif nidus == "spheres" or nidus == "boundingbox" or nidus == "hull":

#             for method in skeletons:
#                 if method == "vmtk" or method == "skeleton":

#                     if polydata is None: #if iḿ looping, keep skeleton

#                         #read img and skeletonise

#                         if method == "vmtk":
#                             img = readNIFTIasVTK(segmfile)
#                             shape = img.GetDimensions()
#                             spacing = img.GetSpacing()
#                             skeoutputpath = os.path.join(ipath, "centerlines", segm, method)
#                             outputpath = os.path.join(ipath, "classification", segm, method, nidus)
#                             network = getvmtkNetwork(img)
#                             polydata = toUniquePointID(network)
                            
#                         else:
#                             img = readNIFTIasSITK(segmfile)
#                             shape = img.GetSize()[::-1]
#                             skeoutputpath = os.path.join(ipath, "centerlines", segm, method, thinmethod)
#                             outputpath = os.path.join(ipath, "classification", segm, method, nidus)
#                             npimgorig = SITKToNumpy(readNIFTIasSITK(os.path.join(ipath, "raw", case+"_0000.nii.gz"))) if grayscale else None
#                             origin = (0,0,0)
#                             spacing = (1,1,1)
#                             img = gaussianSmoothDiscrete(img) #perform gaussian before thinning
#                             npimg = SITKToNumpy(img)
#                             npthinned = binarySegmToBinarySkeleton3(npimg, npimgorig, 0.25)            
#                             ske, _ = binSketoSke2(npthinned, grayscale)
#                             polydata = skeToPolyline(ske)

#                     os.makedirs(skeoutputpath, exist_ok=True)
#                     writeVTKPolydataasVTP(polydata, os.path.join(skeoutputpath, case+".vtp"))

#                     mask = avmMask(shape, origin, spacing, polydata, nidus) 
#                     masksitk = numpyToSITK(mask)
#                     masksitk.SetOrigin(img.GetOrigin())
#                     masksitk.SetSpacing(img.GetSpacing())

#                     os.makedirs(outputpath, exist_ok=True)
#                     writeSITK(masksitk, os.path.join(outputpath, case+".nii.gz"))

#                 else:
#                     print("Invalid skeletonisation method")
#         else:
#             print ("Invalid nidus method")

#     print("--- %s seconds ---" % (time.time() - start_time))



if __name__ == "__main__":
    main()
