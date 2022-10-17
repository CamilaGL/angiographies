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


"""

import vtk
import SimpleITK as sitk
import numpy as np
import argparse
import collections
import math
import time
#from skimage.morphology import convex_hull_image
from angiographies.utils.iositk import readNIFTIasSITK, writeSITK
from angiographies.utils.iojson import writejson, readjson
from angiographies.utils.iovtk import readVTPPolydata
from angiographies.utils.formatconversionsitk import numpyToSITK, SITKToNumpy
from angiographies.utils.imageprocessing import thresholdImageSITK, getLargestConnected




def inside_sphere(_x, _y, _z, center, radius):
    '''Check if a point falls inside a sphere of center and radius'''
    x = np.array([_x, _y, _z]) # The point of interest.
    return (np.linalg.norm(x - center) < radius)

def sphere(mat, rad, c):
    '''Create sphere of rad radius in c center inside mat matrix
    mat: numpy matrix (z y x)
    rad: scalar
    c: tuple with center coordinates (x y z)'''

    #sp = np.zeros((rad*2+1,rad*2+1,rad*2+1))
    spcenter = [rad, rad, rad] #center of the sphere
    for x in range(rad*2+1):
        for y in range(rad*2+1):
            for z in range(rad*2+1):
                #pos = x * incs[0] + y * incs[1] + z * incs[2]
                if inside_sphere(x, y, z, spcenter, rad) and (z+c[2]-rad<mat.shape[0] and y+c[1]-rad<mat.shape[1] and x+c[0]-rad<mat.shape[2]):
                    mat[z+c[2]-rad,y+c[1]-rad,x+c[0]-rad] = 1

def minus(p1, p2):
    return (p1[0]-p2[0],p1[1]-p2[1],p1[2]-p2[2])

def sum(p1, p2):
    return (p1[0]+p2[0],p1[1]+p2[1],p1[2]+p2[2])

def divT(p1, p2): #divide tuple
    return (round(p1[0]/p2[0]),round(p1[1]/p2[1]),round(p1[2]/p2[2]))

def divS(e, p):#scalar to tuple
    return (round(e/p[0]),round(e/p[1]),round(e/p[2]))

def divTS(p, e):#tuple by scalar
    return (p[0]/e,p[1]/e,p[2]/e)

def mulTS(p, e):#tuple by scalar
    return (p[0]*e,p[1]*e,p[2]*e)

def extractNidusSphere(img, radius):
    '''Perform binary closing and opening to extract avm nidus from segmentation
    img: SITKImage
    radius: morphological sphere radius
    returns sitk binary image with nidus voxels True'''
    print("nidus sphere")
    start_time = time.time()
    #img = readNIFTIasSITK(case)
    #imthres = su.getLargestConnected(img) #extract largest island
    npimg = SITKToNumpy(img)
    #print(npimg)
    #binary opening and closing as chenoune2019 does
    vectorRadius = (radius,radius,radius)
    kernel = sitk.sitkBall
    nidusitk = sitk.BinaryMorphologicalClosing(img, vectorRadius, kernel)
    nidusitk = sitk.BinaryMorphologicalOpening(nidusitk, vectorRadius, kernel)
    nidusitk = thresholdImageSITK(nidusitk, 1, 1)
    print(type(npimg[0,0,0]))
    print(type(SITKToNumpy(nidusitk)[0,0,0]))
    #realnidus = np.logical_and(npimg>0,(su.SITKToNumpy(nidusitk))>0)
    realnidus = npimg * SITKToNumpy(nidusitk)
    print(collections.Counter(realnidus[373,152,:])) #255, 0 y 128 para las marcas
    print(collections.Counter(npimg[373,152,:])) #255, 0 y 128 para las marcas
    print(collections.Counter(SITKToNumpy(nidusitk)[373,152,:])) #255, 0 y 128 para las marcas

    #print(realnidus)
    #realnidusnp = realnidus.astype(int)
    #print(realnidusnp)
    realnidusitk = numpyToSITK(realnidus)
    #nidusitk and img to get nidus in image
    #img - (nidusitk and img)
    realnidusitk.SetOrigin(img.GetOrigin())
    realnidusitk.SetSpacing(img.GetSpacing())
    realnidusitk.SetDirection(img.GetDirection())
    #nidusitkconn = getLargestConnected(nidusitk)
    #print(su.SITKToNumpy(nidusitk))
    # extracted = numpyToSITK(np.add(npimg,-realnidus))
    # extracted.SetOrigin(img.GetOrigin())
    # extracted.SetSpacing(img.GetSpacing())
    # extracted.SetDirection(img.GetDirection())
    #writeSITK(extracted,inputf+"AVM_veinarteryCH"+os.path.sep+(((case.split(os.path.sep)[-1]).split(".")[0]))+".nii.gz")
    #writeSITK(realnidusitk,inputf+"AVM_nidusCH"+os.path.sep+(((case.split(os.path.sep)[-1]).split(".")[0]))+".nii.gz")
    print("--- %s seconds ---" % (time.time() - start_time))
    return realnidusitk


#---- change this one to receive centerlines ready

# def extractNidus(inputf, case):
#     print("nidus")
#     img = readNIFTIasSITK(case)
#     #print(img)
#     #connectedFilter = sitk.ConnectedComponentImageFilter()
#     #imgconn = connectedFilter.Execute(img)
#     #imthres = su.thresholdImageSITK(imgconn, 1, 1) #extract largest island
#     imthres = su.getLargestConnected(img) #extract largest island
#     imvtk = NumpyToVTK(SITKToNumpy(imthres)) #sitk to vtk
#     imvtk.SetSpacing(img.GetSpacing()) #vtkimagedata
#     imvtk.SetOrigin(img.GetOrigin())
#     clipped = segmentationToClippedMesh(imvtk,20.0)
#     #binary to mesh and then cut
#     print(clipped)
#     print("done clipping and such")
#     network = centerlines(clipped)
#     print(network)
#     print("done centerlining and such")
#     writeVTKPolydataasVTP(network, inputf+"AVM_centerlines"+os.path.sep+(((case.split(os.path.sep)[-1]).split(".")[0]))+".vtp")

#     # ------- cl and nidus and stuff
#     npimg =SITKToNumpy(imthres)
#     network.BuildLinks()
#     newCl = centerlineediting.removeEndlines(network)
#     newCl2 = centerlinemetrics.getAllMetrics(newCl)
#     print("done removing endlines and such")
#     centerlinemetrics.createTopology(newCl2)
#     print("done creating topology and such")
#     #writeVTKPolydataasVTP(newCl2, inputf+case+"topo.vtp")
#     arrp = centerlineediting.findAVM(newCl2)
#     print("done removing avms and such")
#     #create sphere structuring elements with center and radius as indicated on each arrp entry
#     print(arrp)
#     origin = img.GetOrigin()
#     print(origin)
#     spa = img.GetSpacing()
#     print(spa)
#     mask = np.zeros(npimg.shape)
#     for k in arrp.keys():
#         center = divT(minus(arrp[k][0], origin),spa)
#         rad = divS(arrp[k][1], spa)
#         sphere(mask, rad[0], center)
#     print(type(mask))
#     print(type(npimg))
#     nidus = np.logical_and(npimg>0, mask>0)
#     #print(nidus.astype(int))
#     nidusitk = su.numpyToSITK(nidus.astype(int))
#     nidusitk.SetOrigin(img.GetOrigin())
#     nidusitk.SetSpacing(img.GetSpacing())
#     nidusitk.SetDirection(img.GetDirection())
#     nidusitkconn = su.getLargestConnected(nidusitk)
#     extracted = su.numpyToSITK(np.add(npimg,-(su.SITKToNumpy(nidusitkconn))))
#     extracted.SetOrigin(img.GetOrigin())
#     extracted.SetSpacing(img.GetSpacing())
#     extracted.SetDirection(img.GetDirection())
#     io.writeSITK(extracted,inputf+"AVM_veinarteryBA"+os.path.sep+(((case.split(os.path.sep)[-1]).split(".")[0]))+".nii.gz")
#     io.writeSITK(nidusitkconn,inputf+"AVM_nidusBA"+os.path.sep+(((case.split(os.path.sep)[-1]).split(".")[0]))+".nii.gz")


def getAdjacentPoints(skeleton, v):
    '''Creates a list with the points connected to v
    skeleton: polydata with lines and points
    v: point id'''
    adjacent = []
    cellIds = vtk.vtkIdList()
    skeleton.GetPointCells(v, cellIds) #get cells incident
    incident = [cellIds.GetId(c) for c in range(cellIds.GetNumberOfIds())]
    for c in incident:
        if (skeleton.GetPoint(skeleton.GetCell(c).GetPointId(0)) == skeleton.GetPoint(v)): #first point, find next point
            adjacent.append(skeleton.GetPoint(skeleton.GetCell(c).GetPointId(1)))
        elif (skeleton.GetPoint(skeleton.GetCell(c).GetPointId(skeleton.GetCell(c).GetNumberOfPoints()-1)) == skeleton.GetPoint(v)): #last point, find second to last point
            adjacent.append(skeleton.GetPoint(skeleton.GetCell(c).GetPointId(skeleton.GetCell(c).GetNumberOfPoints()-2)))
    return adjacent


def getAdjacentPointsExtended(skeleton, v, arrayname = "Radius"):
    '''Creates a list with the points connected to v, but considering the line extended according to radius | edt
    skeleton: polydata with lines and points
    v: point id'''
    adjacent = []
    cellIds = vtk.vtkIdList()
    p1 = skeleton.GetPoint(v)
    skeleton.GetPointCells(v, cellIds) #get cells incident
    incident = [cellIds.GetId(c) for c in range(cellIds.GetNumberOfIds())]
    for c in incident:
        if (skeleton.GetPoint(skeleton.GetCell(c).GetPointId(0)) == skeleton.GetPoint(v)): #first point, find next point
            p2id = skeleton.GetCell(c).GetPointId(1)
        elif (skeleton.GetPoint(skeleton.GetCell(c).GetPointId(skeleton.GetCell(c).GetNumberOfPoints()-1)) == skeleton.GetPoint(v)): #last point, find second to last point
            p2id = skeleton.GetCell(c).GetPointId(skeleton.GetCell(c).GetNumberOfPoints()-2)
        p2 = skeleton.GetPoint(p2id)
        lenAB = math.sqrt(pow(p1[0] - p2[0], 2.0) + pow(p1[1] - p2[1], 2.0) + pow(p1[2] - p2[2], 2.0))
        vectAB = minus(p2, p1)
        normAB = divTS(vectAB, lenAB)
        extendto = skeleton.GetPointData().GetArray(arrayname).GetValue(p2id)
        p3 = sum(mulTS(normAB, extendto), p1)
        adjacent.append(p3)
    return adjacent


def getFurthestDist(skeleton, v):
    '''skeleton: polydata with lines and points
    v: point id'''
    furthestDist = 0
    cellIds = vtk.vtkIdList()
    skeleton.GetPointCells(v, cellIds) #get cells incident
    incident = [cellIds.GetId(c) for c in range(cellIds.GetNumberOfIds())]
    for c in incident:
        if (skeleton.GetPoint(skeleton.GetCell(c).GetPointId(0)) == skeleton.GetPoint(v)): #first point, find distance to next point
            dist = math.sqrt(vtk.vtkMath.Distance2BetweenPoints(skeleton.GetPoint(skeleton.GetCell(c).GetPointId(1)), skeleton.GetPoint(skeleton.GetCell(c).GetPointId(0))))
            if dist > furthestDist:
                furthestDist = dist
        elif (skeleton.GetPoint(skeleton.GetCell(c).GetPointId(skeleton.GetCell(c).GetNumberOfPoints()-1)) == skeleton.GetPoint(v)):
            dist = math.sqrt(vtk.vtkMath.Distance2BetweenPoints(skeleton.GetPoint(skeleton.GetCell(c).GetPointId(skeleton.GetCell(c).GetNumberOfPoints()-2)), skeleton.GetPoint(skeleton.GetCell(c).GetPointId(skeleton.GetCell(c).GetNumberOfPoints()-1))))
            if dist > furthestDist:
                furthestDist = dist

    return furthestDist

def areAdjacent(skeleton, v1, v2):
    cellIds1 = vtk.vtkIdList()
    cellIds2 = vtk.vtkIdList()
    skeleton.GetPointCells(v1, cellIds1) #get cells incident
    skeleton.GetPointCells(v2, cellIds2) #get cells incident
    v1cells = [cellIds1.GetId(i) for i in range(cellIds1.GetNumberOfIds())]
    commonCells = [cellIds2.GetId(i) for i in range(cellIds2.GetNumberOfIds()) if cellIds2.GetId(i) in v1cells]
    if len(commonCells)>0:
        return True
    return False

def multipleSpheres(imgshape, spheres, origin=(0,0,0), spacing=(1,1,1)):
    '''imgshape: Image shape to create mask.
    spheres: dictionary with tuples where the first element is sphere center and second is sphere radius.
    origin: if sphere coords are in real world, indicate origin to translate to voxels.
    spacing: if sphere coords are in real world, indicate voxel spacing to translate to voxels.'''
    print(imgshape)
    print(origin)
    origin=(0,0,0) #origin always 0, because skeleton starts at 0 and vmtk seems to be ignoring the origin value?
    mask = np.zeros(imgshape)
    for k in spheres.keys():
        print(k)
        center = divT(minus(spheres[k][0], origin),spacing)
        print(spheres[k][0])
        print(center)
        rad = divS(spheres[k][1], spacing)
        print(rad)
        sphere(mask, rad[0], center)
    return mask

def convexHull(imgshape, coords, origin=(0,0,0), spacing=(1,1,1)):
    '''image is a binary np matrix
    returns a binary image with true on the voxels that mark the center of the spiders and their adjacent points'''
    print("convex hull delimiters")
    mask = np.zeros(imgshape)
    for k in coords:
        for p in coords[k]:
            mask[(divT(minus(p, origin),spacing))[::-1]] = 1
    # hull = convex_hull_image(mask)

    return mask


def findAVMSpiders(skeleton, method="spheres"):
    '''skeleton: vtkpolydata with unique points (each location has its own id)
    method: how to report the spiders (spheres: max radius for each sphere center, hull: all adjacent points to get convex hull)'''

    #create dictionary with points that could potentially belong to the avm and count references
    vertexDict = {}
    maxdeg=0
    avm=0
    for i in range(skeleton.GetNumberOfPoints()):  
        cellIds = vtk.vtkIdList()
        skeleton.GetPointCells(i, cellIds) #get cells adjacent to our point
        degree = cellIds.GetNumberOfIds() #how many cells are incident to this point?
        if degree > 5: #We consider this a potential avm point
            vertexDict[i] = degree  #I save this one
        if degree > maxdeg:
            avm=i #This one has the highest degree until now
            maxdeg=degree
    print(vertexDict)    
    #get potential avm nodes adjacent to the main node (highest deg) and their point coords
    avmPoints = {}
    avmPoints[avm] = [skeleton.GetPoint(avm)]#.append(avm)
    for k in vertexDict:
        if areAdjacent(skeleton, avm, k) and k not in avmPoints: #check we haven't added it yet bc we'd have twice the same --- Not necessary anymore
            #avmNodes.append(k)
            avmPoints[k] = [skeleton.GetPoint(k)]

    #find max radius for sphere centered on each node
    if method=="spheres":
        for k in avmPoints:
            avmPoints[k].append(getFurthestDist(skeleton, k))
    elif method == "hull":
        for k in avmPoints:
            avmPoints[k].extend(getAdjacentPoints(skeleton, k))
        #get all points connected to sphere centers which will allow to find convex hull later
    #return central points and sphere radius

    return avmPoints
    

def avmMask(imgshape, origin, spacing, skeleton, method="spheres"):
    '''skeleton: vtkpolydata with unique points (each location has its own id)
    method: how to report the spiders (spheres: max radius for each sphere center, hull: all adjacent points to get convex hull)'''
    print("doing spiders with method", method)
    spiders = findAVMSpiders(skeleton, method)
    print(spiders)
    print(method)
    if method=="spheres":
        print("doing spheres")
        mask = multipleSpheres(imgshape, spiders, origin, spacing)
    elif method == "hull":
        mask = convexHull(imgshape, spiders, origin, spacing)
    return mask


# def processCenterline(inputf, case):
#     img = io.readNIFTIasSITK(inputf+"AVM_whole"+os.path.sep+(((case.split(os.path.sep)[-1]).split(".")[0]))+".nii.gz")
#     print(img.GetSize())
#     imgisland = su.getLargestConnected(img)
#     npimg =su.SITKToNumpy(imgisland)
#     network = io.readVTPPolydata(case)
#     network.BuildLinks()
#     newCl = centerlineediting.removeEndlines(network)
#     newCl2 = centerlinemetrics.getAllMetrics(newCl)
#     print("done removing endlines and such")
#     centerlinemetrics.createTopology(newCl2)
#     print("done creating topology and such")
#     #writeVTKPolydataasVTP(newCl2, inputf+case+"topo.vtp")
#     arrp = centerlineediting.findAVM(newCl2)
#     print("done removing avms and such")
#     #create sphere structuring elements with center and radius as indicated on each arrp entry
#     print(arrp)
#     origin = img.GetOrigin()
#     print(origin)
#     spa = img.GetSpacing()
#     print(spa)
#     mask = np.zeros(npimg.shape)
#     for k in arrp.keys():
#         print(k)
#         center = divT(minus(arrp[k][0], origin),spa)
#         rad = divS(arrp[k][1], spa)
#         sphere(mask, rad[0], center)
#     print(type(mask))
#     print(type(npimg))
#     nidus = np.logical_and(npimg>0, mask>0)
#     #print(nidus.astype(int))
#     nidusitk = su.numpyToSITK(nidus.astype(int))
#     nidusitk.SetOrigin(img.GetOrigin())
#     nidusitk.SetSpacing(img.GetSpacing())
#     nidusitk.SetDirection(img.GetDirection())
#     nidusitkconn = su.getLargestConnected(nidusitk)
#     extracted = su.numpyToSITK(np.add(npimg,-(su.SITKToNumpy(nidusitkconn))))
#     extracted.SetOrigin(img.GetOrigin())
#     extracted.SetSpacing(img.GetSpacing())
#     extracted.SetDirection(img.GetDirection())
#     io.writeSITK(extracted,inputf+"AVM_veinartery"+os.path.sep+(((case.split(os.path.sep)[-1]).split(".")[0]))+".nii.gz")
#     io.writeSITK(nidusitkconn,inputf+"AVM_nidus"+os.path.sep+(((case.split(os.path.sep)[-1]).split(".")[0]))+".nii.gz")
        

def main():

    parser = argparse.ArgumentParser()
    #parser.add_argument("-s", help="seed acquisition method", default="5", required=False)
    parser.add_argument("-ifile", help="path to segmentation or centerline", default="", required=True)
    parser.add_argument("-ofile", help="path to output mask", default="", required=True)
    parser.add_argument("--sphere", help="process volume with chenoune sphere", action="store_true", required=False, default=False)
    parser.add_argument("-rad", help="sphere radius for morphological extraction", type=int, default=30, required=False)
    parser.add_argument("--spidersphere", help="process volume with babin skeleton", action="store_true", required=False, default=False)
    parser.add_argument("--spiderhull", help="process volume with babin skeleton", action="store_true", required=False, default=False)


    args = parser.parse_args()
    inputf = args.ifile
    outputf = args.ofile
    sphere = args.sphere
    spidersphere = args.spidersphere
    spiderhull = args.spiderhull
    sphereradius = args.rad

    if sphere:
        print("morphological spheres")
        img = readNIFTIasSITK(inputf)
        print(img.GetOrigin())
        print(img.GetSize())
        print(SITKToNumpy(img).shape)
        nidus = extractNidusSphere(img, sphereradius)#spheres
        writeSITK(nidus, outputf)
    else:
        img = readVTPPolydata(inputf)
        #print(inputf[:-19]+".json")
        infodict = readjson(inputf[:-19]+".json") if "vmtk" in inputf else readjson(inputf[:-12]+".json")
        print(infodict)
        method = "spheres" if spidersphere else "hull"
        numpyshape = infodict["shape"][::-1]
        numpyorigin = (0,0,0)
        numpyspa = infodict["spacing"] if "vmtk" in inputf else (1,1,1)

        mask = avmMask(numpyshape, numpyorigin, numpyspa, img, method)
        masksitk = numpyToSITK(mask)
        masksitk.SetOrigin(infodict["origin"])
        masksitk.SetSpacing(infodict["spacing"])
        print(masksitk.GetOrigin())
        print(masksitk.GetSize())
        writeSITK(masksitk, outputf)



if __name__ == "__main__":
    main()
