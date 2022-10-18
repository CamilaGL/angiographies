"""
This is based on the thinning algorithm described here 10.1109/DICTA.2009.13 MINUS the templates (potential TO DO?)
And the voxel consideration order is given by the squared euclidean distance to the nearest background voxel
Exceptionally not optimal

Running environment requirements: 

    numpy
    sitk

"""

import numpy as np
import argparse
import time
#from scipy import ndimage
from angiographies.utils.iositk import readNIFTIasSITK, writeSITK
from angiographies.utils.formatconversionsitk import numpyToSITK, SITKToNumpy
from angiographies.utils.imageprocessing import binClosing

# --------- ordered thinning with different criterions ---------------

def binarySegmToBinarySkeleton(img, inputgrayscale=None):
    '''Binary thinning. inputgrayscale is required if thinning order is grayvalue. Else, order is euclidean distance to background.
    img: sitk image
    inputgrayscale: sitk image
    returns sitk image'''
    
    npimg = SITKToNumpy(img)
    #npimg = np.array([[[0,1,1],[0,1,1],[1,1,1]],[[0,1,0],[0,1,1],[0,1,1]],[[0,1,0],[0,1,0],[1,1,1]]], dtype=np.ubyte)
    weighted = np.zeros(npimg.shape, dtype=np.intc) #here we're going to be assigning the order priority to our foreground voxels
   
    #weighted = None
    start_time = time.time()
    if inputgrayscale is None:
        print("Euclidean weighting")
        #weighted = np.empty(npimg.shape, dtype=np.intc) #here we're going to be assigning the order priority to our foreground voxels
        #print(np.iinfo(np.intc).max)
        #weighted[:,:,:] = np.iinfo(np.intc).max
        #weighted = np.zeros(npimg.shape, dtype=np.intc) #here we're going to be assigning the order priority to our foreground voxels
        getWeightedImageEuclidean(npimg, weighted) #we're editing weighted here!
    else:
        print("Grayscale weighting")
        npimgorig = SITKToNumpy(inputgrayscale)
        #weighted = np.empty(npimg.shape, dtype=np.intc) #here we're going to be assigning the order priority to our foreground voxels
        #print(np.iinfo(np.intc).max)
        #weighted[:,:,:] = np.iinfo(np.intc).max
        if npimg.shape == npimgorig.shape: #check segmentation and grayscale images have same dimentions
            getWeightedImageGrayscale(npimg, npimgorig, weighted) #we're editing weighted here!
    print("Weighted image ready.")
    newOrderedThinning(npimg, weighted) #we're editing npimg here!
    print("--- %s seconds ---" % (time.time() - start_time))
    
    thinnedsitk = numpyToSITK(npimg)
    thinnedsitk.SetOrigin(img.GetOrigin())
    thinnedsitk.SetSpacing(img.GetSpacing())
    thinnedsitk.SetDirection(img.GetDirection())
    return thinnedsitk


def intersection(lst1, lst2):
    '''Returns the intersection between two lists'''
    temp = set(lst2)
    return [value for value in lst1 if value in temp]

def getNNeighIndexes(v, maxlim, n=26, includev=False):
    '''Get a list with all the indexes of v's 6/18/26 connected neighbourhood unrolled.
    v is the voxel to which we're computing the neighbourhood (in tuple format).
    maxlim is the shape of the volume.
    n is the number of neigbhours we're computing.
    includev indicates if we want the voxel of interest included in out comprehensive neighbourhood list'''
    list_max = [x + y for x, y in zip(v, [1,1,1])]
    list_min = [x - y for x, y in zip(v, [1,1,1])]
    neigh_min = [max((line_min, 0)) for line_min in list_min] #make sure the min index is still inside the volume
    neigh_max = [min((line_max, line_shape-1)) for line_max, line_shape in zip(list_max, maxlim)] #make sure the max index is still inside the volume
    #now unroll
    indexes = []
    for i in range(neigh_min[0], neigh_max[0]+1): 
        for j in range(neigh_min[1], neigh_max[1]+1):
            for k in range(neigh_min[2], neigh_max[2]+1):
                if n == 26:
                    indexes.append((i,j,k)) #all are valid
                elif n == 18:
                    if (abs(v[0]-i) + abs(v[1]-j) + abs(v[2]-k)) <=2: #can only move in two directions at a time
                        indexes.append((i,j,k))
                elif n == 6:
                    if (abs(v[0]-i) + abs(v[1]-j) + abs(v[2]-k)) <= 1: #can only move in one direction at a time
                        indexes.append((i,j,k))
    if not includev:
        indexes.remove(v) #I don't want the voxel of interest to be included
    return indexes

def isBoundary(img, v):
    '''Check if the 26-neighbourhood of v has at least on background voxel'''
    indexes26 = getNNeighIndexes(v, img.shape, 26, False)
    boundary = False
    for p in indexes26[:]:
        if img[p] == 0:
            boundary = True #found one background
    return boundary

def isBackgroundSimple(img, v):
    '''Check if the voxels of *value* in the *neigh* neighbourhood of *v* voxel are a connected component'''
    indexes18 = getNNeighIndexes(v, img.shape, 18, False) #get neighbours without v
    indexes6 = getNNeighIndexes(v, img.shape, 6, False)
    #nz = list(zip(*(np.where(img==value)))) #get img where equals value
    for p in indexes18[:]:
        if img[p] != 0:
            indexes18.remove(p) #remove from the index list all neighbours that are not background
    
    if len(intersection(indexes18, indexes6))==0: #there are no background in the 6 connected
        return False
    
    else: #checked that part, now check if 18-neigh background is connected (flood fill)
        visit = [indexes18.pop()]
        while len(visit) != 0: #go over all the connected neighbours
            p = visit.pop()
            #print(p)
            for e in indexes18[:]: #check if any of the not visited is a neighbour (iterate over a copy so we can delete elements)
                #print(e)
                if max(abs(p[0]-e[0]),abs(p[1]-e[1]),abs(p[2]-e[2])) == 1:
                    visit.append(e)
                    indexes18.remove(e)
        if len(indexes18) == 0: #I visited all the background voxels in one connected pass
            return True
        else:
            return False

def isForegroundSimple(img, v):
    '''Check if the voxels of foreground in the 26 neighbourhood of *v* voxel are a connected component
    (which means, the removal of v does not affect the foreground connectivity)'''
    indexes = getNNeighIndexes(v, img.shape, 26, False) #get 26 connected component without v
    numneigh = len(indexes)  
    for p in indexes[:]:
        if img[p] != 1: #remove from the index list all neighbours that are not foreground
            indexes.remove(p)
    
    #if no foreground, v is not simple. If only one foreground, v is endline. If all neigh are foreground, v is not simple.
    if len(indexes)<=1 or len(indexes) == numneigh: 
        return False
    
    visit = [indexes.pop()] #pick one voxel to study connectivity
    while len(visit) != 0: #go over all the connected neighbours (flood fill)
        p = visit.pop()
        #print(p)
        for e in indexes[:]: #check if any of the not visited is a neighbour (iterate over a copy so we can delete elements in indexes)
            #print(e)
            if max(abs(p[0]-e[0]),abs(p[1]-e[1]),abs(p[2]-e[2])) == 1:
                visit.append(e)
                indexes.remove(e)
    if len(indexes) == 0: #I visited all the foreground voxels in one connected pass
        return True
    else:
        return False


def isConnectedComponent(img, value=1):
    '''Check if the elements of *value* inside the volume *img* are a connected component'''
    nz = list(zip(*(np.where(img==value)))) #get img where equals value
    if len(nz)==0:
        return False
    #print(nz)
    visit = [nz.pop()]
    while len(visit) != 0: #go over all the connected neighbours
        p = visit.pop()
        #print(p)
        for e in nz[:]: #check if any of the not visited is a neighbour (iterate over a copy so we can delete elements in nz)
            #print(e)
            if max(abs(p[0]-e[0]),abs(p[1]-e[1]),abs(p[2]-e[2])) == 1:
                visit.append(e)
                nz.remove(e)
    if len(nz) == 0: #I visited all the *value* voxels in one connected pass
        return True
    else:
        return False

def isRedundant(img, v):
    if isForegroundSimple(img, v) and isBackgroundSimple(img, v):
        return True
    else:
        return False

# def orderedThinning(img, weightedImg):
#     nz = img != 0
#     #print(nz)
#     cleanednz = np.where(nz, weightedImg, np.nan)
#     #print(cleanednz)
#     orderednz = np.argsort(cleanednz, axis=None)[:nz.sum()]
#     #print(orderednz)
#     ind = np.unravel_index(orderednz, img.shape)
#     #print(ind)
#     #print(weightedImg[ind])
#     #print(img)
#     for x,y,z in zip(*ind): #iterate over ordered nonzero voxels
#         #print("mirando")
#         if isRedundant(img, (x,y,z)):
#             #print(img)
#             #print("este ", x, y, z, "es redundante")
#             img[x,y,z] = 0
#     #print(img)

def newOrderedThinning(img, weightedImg):
    '''Ordered thinning algorithm considering the weightedImg values as order priority.
    Loop implemented according to the She et al (2009) paper. (Templates NOT implemented)'''
    nz = img != 0
    #print(nz)
    cleanednz = np.where(nz, weightedImg, np.nan)
    #print(cleanednz)
    orderednz = np.argsort(cleanednz, axis=None)[:nz.sum()]
    #print(orderednz)
    ind = np.unravel_index(orderednz, img.shape)
    #print(ind)
    #print(weightedImg[ind])
    #print(img)
    notmarked = list(zip(*ind)) #get all foreground voxels
    marked = []
    stillworking = True
    while stillworking:
        stillworking=False
        for t in notmarked[:]: #iterate over ordered nonzero voxels
            #print("mirando")
            if isBoundary(img, t):
                marked.append(t)
                notmarked.remove(t)
        stillmarked = True
        while stillmarked:
            stillmarked = False
            for t in marked[:]:
                if isRedundant(img, t):
                    img[t] = 0
                    marked.remove(t)
                    stillmarked=True #did something, let's keep looping
                    stillworking=True #did something, let's keep looping
            if len(marked)==0:
                stillmarked=False #no more marked to check
        for t in marked[::-1]:
            notmarked.insert(0,t)
        marked = []
    #print(img)


# def getWeightedImageEuclideanScipy(img):
#     '''Generate a new image where the foreground voxels are weighted according to the squared euclidean distance
#     (relevant for their ordered consideration to perform a binary thinning)'''
#     profiled = ndimage.distance_transform_edt(img)
#     return profiled


def getWeightedImageEuclidean(img, profiled):
    '''Generate a new image where the foreground voxels are weighted according to the squared euclidean distance
    (relevant for their ordered consideration to perform a binary thinning)'''
    specified = np.nonzero(img) #get nonzero voxels
    for (x,y,z) in zip(*specified): #iterate over nonzero voxels
        profiled[x,y,z] = getMinEuclideanToBackground(img, (x,y,z))
        #profiled[x,y,z] = getProfileMeasure(img, (x,y,z), getMinEuclideanToBackground(img, (x,y,z)))

def getWeightedImageGrayscale(img, original, profiled):
    '''Generate a new image where the foreground voxels are assigned their (original) grayscale value
    (relevant for their ordered consideration to perform a binary thinning)'''
    specified = np.nonzero(img) #get nonzero voxels in segmentation
    for (x,y,z) in zip(*specified): #iterate over nonzero voxels
        profiled[x,y,z] = original[x,y,z]

def getWeightedProfileMeasure(img, profiled):
    '''Generate a new image where the foreground voxels are weighted according to the profile measure proposed by BABIN Thesis
    (relevant for their ordered consideration to perform a binary thinning)'''
    specified = np.nonzero(img) #get nonzero voxels
    for (x,y,z) in zip(*specified): #iterate over nonzero voxels
        profiled[x,y,z] = getProfileMeasure(img, (x,y,z), getMinEuclideanToBackground(img, (x,y,z)))

def inside_sphere(_x, _y, _z, center, radius):
    '''Check if x,y,z is inside the sphere of *center* and *radius*'''
    x = np.array([_x, _y, _z]) # The point of interest.
    return (np.linalg.norm(x - center) <= radius)

def getProfileMeasure(img, v, delta):
    '''Count the number of foreground voxels included in a sphere of radius square root *delta* and center in v (excluding v)'''
    sqrtdelta = np.sqrt(delta)
    #print(sqrtdelta)
    intsqrtdelta = int(np.ceil(sqrtdelta))
    #print(intsqrtdelta)
    neigh_min, neigh_max = getNeighLimits(img, v, [intsqrtdelta, intsqrtdelta, intsqrtdelta]) #get cube of side 2*ceil(sqrt(delta))
    #print ("limits are  min", neigh_min, "and max", neigh_max)
    profilevalue = -1 #instead of zero, to account for the center (voxel of interest) that we have to ignore
    neighbours = np.nonzero(img[neigh_min[0]:neigh_max[0]+1,neigh_min[1]:neigh_max[1]+1,neigh_min[2]:neigh_max[2]+1])
    # for i in range(neigh_min[0], neigh_max[0]+1):
    #     for j in range(neigh_min[1], neigh_max[1]+1):
    #         for k in range(neigh_min[2], neigh_max[2]+1):
    for i,j,k in zip(*neighbours):
        #print("neighbour is ", i,j,k)
        if inside_sphere(i+neigh_min[0],j+neigh_min[1],k+neigh_min[2],v,sqrtdelta):# and img[i,j,k]==1: #distance to voxel is inside sphere (cube is generated to make stuff easier)
            profilevalue = profilevalue + 1
            #print("is inside")
    return profilevalue

def getNeighLimits(img, v, radius): #v as tuple, radius as list
    '''For certain radius value, get valid indexes in img of a cube with center in v'''
    list_max = [x + y for x, y in zip(v, radius)]
    list_min = [x - y for x, y in zip(v, radius)]
    neigh_min = [max((line_min, 0)) for line_min in list_min]
    neigh_max = [min((line_max, line_shape-1)) for line_max, line_shape in zip(list_max, img.shape)]
    return neigh_min, neigh_max

def stretchMask(img, v, radius): #v as tuple, radius as list
    '''Stretch neighbourhood mask radius and get valid limits inside img with center in v'''
    newradius = [x + 1 for x in radius] 
    neigh_min, neigh_max = getNeighLimits(img, v, newradius)
    return neigh_min, neigh_max, newradius

def getMinEuclideanToBackground(img, v):
    '''Returns the minimal squared Euclidean distance (delta) from foreground voxel v to the background (nearest background voxel)'''
    radius = [0,0,0]
    found = False
    delta = np.dot(np.asarray(img.shape).T,np.asarray(img.shape)) #distance from first to last voxel (max distance between two voxels in this volume)
    while not found: #look for a background voxel
        neigh_min, neigh_max, radius = stretchMask(img, v, radius) #look on a bigger neighbourhood
        neighbourhood = img[neigh_min[0]:neigh_max[0]+1,neigh_min[1]:neigh_max[1]+1,neigh_min[2]:neigh_max[2]+1] #get submatrix corresponding to neighbourhood
        background = np.where(neighbourhood==0) #get background surrounding our voxel
        
        if len(list(zip(*background)))>0: #I have found background
            
            found = True
            for (x,y,z) in zip(*background):
                vector = np.asarray([vox-(back+minbias) for vox, back, minbias in zip(v,[x,y,z],neigh_min)]) #get distance vector
                
                newdelta = np.dot(vector.T, vector) #Euclidean distance squared
                if newdelta < delta: delta = newdelta #keep min distance, I don't need to keep the voxel indexes
            
    return delta


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-ifile", help="path to case segmentation", default="", required=True)
    parser.add_argument("-ofile", help="path to output folder+case", default="", required=True)
    parser.add_argument("-gsfile", help="path to original grayscale image", default=None, required=False)
    parser.add_argument("--closing", help="perform morphological closing", action="store_true", required=False, default=False)
    parser.add_argument("-rp", help="repeat closing", default=1, type=int, required=False)
    parser.add_argument("-rad", help="repeat closing", default=1, type=int, required=False)
    #parser.add_argument("--connected", help="extract connected components", action="store_true", required=False, default=False)


    args = parser.parse_args()
    inputf = args.ifile
    outputf = args.ofile
    inputgrayscale = args.gsfile
    closing = args.closing
    rp = args.rp
    rad = args.rad
    #conncomp = args.connected
    print("Starting thinning")
    img = readNIFTIasSITK(inputf)
    imgorig = readNIFTIasSITK(inputgrayscale) if inputgrayscale is not None else None
    if closing:
        for _ in range(rp):
            img = binClosing(img, 2)
    thinnedsitk = binarySegmToBinarySkeleton(img, imgorig)

    writeSITK(thinnedsitk,outputf)



if __name__ == "__main__":
    main()

