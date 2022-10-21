"""
Based on the thinning algorithm loop described by Larrabide 2007 
Voxel consideration order is given by the squared euclidean distance to the nearest background voxel
edt: https://github.com/seung-lab/euclidean-distance-transform-3d

Running environment requirements: 

    numpy
    sitk
    edt

"""


import numpy as np
import argparse
import time
import edt
import heapq
from angiographies.utils.iositk import readNIFTIasSITK, writeSITK
from angiographies.utils.formatconversionsitk import numpyToSITK, SITKToNumpy
from angiographies.utils.imageprocessing import binClosing, gaussianSmoothDiscrete

# --------- ordered thinning with different criterions ---------------

def binarySegmToBinarySkeleton(img):
    '''Binary thinning. Order is squared euclidean distance to background.
    img: sitk image
    returns sitk image'''
    
    npimg = SITKToNumpy(img)
    #npimg = np.array([[[0,1,1],[0,1,1],[1,1,1]],[[0,1,0],[0,1,1],[0,1,1]],[[0,1,0],[0,1,0],[1,1,1]]], dtype=np.ubyte)
    #weighted = np.zeros(npimg.shape, dtype=np.intc) #here we're going to be assigning the order priority to our foreground voxels
    #weighted = None
    start_time = time.time()
    print("Euclidean weighting")
    weighted = getWeightedImageEuclideanDistanceTransform(npimg)
    print("--- %s seconds ---" % (time.time() - start_time))

    newOrderedThinning2(npimg, weighted) #we're editing npimg here!
    print("--- %s seconds ---" % (time.time() - start_time))
    
    thinnedsitk = numpyToSITK(npimg)
    thinnedsitk.SetOrigin(img.GetOrigin())
    thinnedsitk.SetSpacing(img.GetSpacing())
    thinnedsitk.SetDirection(img.GetDirection())
    return thinnedsitk


def binarySegmToBinarySkeleton3(npimg, npimgorig = None, thresh=None):
    '''Binary thinning. Order is squared euclidean distance to background.
    img: numpy image
    returns numpy image'''
    #print(npimg.dtype)
    if thresh is not None:
        npimg[npimg<=thresh]=0
        npimg[npimg>thresh]=1
        npimg = npimg.astype(np.int8)
    #npimg = np.array([[[0,1,1],[0,1,1],[1,1,1]],[[0,1,0],[0,1,1],[0,1,1]],[[0,1,0],[0,1,0],[1,1,1]]], dtype=np.ubyte)
    #weighted = np.zeros(npimg.shape, dtype=np.intc) #here we're going to be assigning the order priority to our foreground voxels
    #npimgpadded = np.pad(npimg, 1)
    #print("orig size", npimg.shape)
    #print("padded size", npimgpadded.shape)
    #weighted = None
    start_time = time.time()
    npimgpadded = np.pad(npimg, 1)
    if npimgorig is None:
        print("Euclidean weighting")
        weighted = getWeightedImageEuclideanDistanceTransform(npimgpadded)
        print("--- %s seconds ---" % (time.time() - start_time))
        #print("weighted size", weighted.shape)
    else:
        print("Grayscale weighting")
        minval = np.amin(npimgorig)
        print(npimgorig[1,1,1])
        print(minval)
        npimgorig = npimgorig - minval + 1
        print(npimgorig[1,1,1])
        npimgorigpadded = np.pad(npimgorig, 1)
        #weighted = np.pad(weighted, 1)
        if npimg.shape == npimgorig.shape: #check segmentation and grayscale images have same dimentions
            weighted = np.where(npimgpadded, npimgorigpadded, 0)
        print("--- %s seconds ---" % (time.time() - start_time))
    newOrderedThinning6(npimgpadded, weighted) #we're editing npimg here!
    print("--- %s seconds ---" % (time.time() - start_time))
    
    distske = np.where(npimgpadded, weighted, 0)

    return distske[1:-1,1:-1,1:-1]


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


def getNNeighIndexes2(v, n=26, includev=False):
    '''Get a list with all the indexes of v's 6/18/26 connected neighbourhood unrolled. We assume a padded image and don't check out of bounds.
    v is the voxel to which we're computing the neighbourhood (in tuple format).
    n is the number of neigbhours we're computing.
    includev indicates if we want the voxel of interest included in out comprehensive neighbourhood list'''

    indexes = []
    for i in range(-1,2): 
        for j in range(-1,2):
            for k in range(-1,2):
                if n == 26:
                    indexes.append((v[0]+i,v[1]+j,v[2]+k)) #all are valid
                elif n == 18:
                    if (abs(i) + abs(j) + abs(k)) <=2: #can only move in two directions at a time
                        indexes.append((v[0]+i,v[1]+j,v[2]+k))
                elif n == 6:
                    if (abs(i) + abs(j) + abs(k)) <= 1: #can only move in one direction at a time
                        indexes.append((v[0]+i,v[1]+j,v[2]+k))
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


def isBoundary2(img, v):
    '''Check if the 26-neighbourhood of v has at least on background voxel'''
    for i in range(-1, 2): 
        for j in range(-1, 2):
            for k in range(-1, 2):
                if img[v[0]+i,v[1]+j,v[2]+k] == 0:
                    return True
    return False


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


def isBackgroundSimple2(img, v):
    '''Check if the voxels of *value* in the *neigh* neighbourhood of *v* voxel are a connected component'''
    indexes18 = getNNeighIndexes2(v, 18, False) #get neighbours without v
    indexes6 = getNNeighIndexes2(v, 6, False)
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


def isForegroundSimple2(img, v):
    '''Check if the voxels of foreground in the 26 neighbourhood of *v* voxel are a connected component
    (which means, the removal of v does not affect the foreground connectivity)'''
    numneigh=0
    indexes=[]
    for i in range(-1, 2): 
        for j in range(-1, 2):
            for k in range(-1, 2):
                if img[v[0]+i,v[1]+j,v[2]+k] == 1:
                    numneigh=numneigh+1
                    if (abs(i) + abs(j) + abs(k)) != 0: #we're excluding the center point v
                        indexes.append((v[0]+i,v[1]+j,v[2]+k))

    #if no foreground, v is not simple. If only one additional foreground, v is endline. If all neigh are foreground, v is not simple.
    if numneigh<=2 or numneigh == 27: 
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

def isEndPoint2(img, v):
    '''Check if the voxels of foreground in the 26 neighbourhood of *v* voxel are a connected component
    (which means, the removal of v does not affect the foreground connectivity)'''
    numneigh=0
    for i in range(-1, 2): 
        for j in range(-1, 2):
            for k in range(-1, 2):
                if img[v[0]+i,v[1]+j,v[2]+k] == 1:
                    numneigh=numneigh+1

    #if no foreground, v is not simple. If only one additional foreground, v is endline. If all neigh are foreground, v is not simple.
    if numneigh<=2: 
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

def isRedundant2(img, v):
    if isForegroundSimple2(img, v) and isBackgroundSimple2(img, v):
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


def newOrderedThinning2(img, weightedImg):
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
    print("Now looping")
    start_time = time.time()
    while stillworking:
        outerloop_start_time = time.time()
        stillworking=False
        for t in notmarked[:]: #iterate over ordered nonzero voxels
            #print("mirando")
            if isBoundary(img, t):
                if isRedundant(img,t):
                    marked.append(t)
                    notmarked.remove(t)
        print("--- selecting redundant boundary %s seconds ---" % (time.time() - outerloop_start_time))
        stillmarked = True
        innerloop_start_time = time.time()
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
        print("--- removing redundant boundary %s seconds ---" % (time.time() - innerloop_start_time))
        for t in marked[::-1]:
            notmarked.insert(0,t)
        marked = []
        print("--- one whole pass %s seconds ---" % (time.time() - outerloop_start_time))
    #print(img)
    print("--- the complete loop took %s seconds ---" % (time.time() - start_time))


def newOrderedThinning3(img, weightedImg):
    '''Ordered thinning algorithm considering the weightedImg values as order priority.
    Loop implemented according to the Larrabide 2007 paper.
    Optimised with heapq and some function calling reduction'''
    nz = np.nonzero(img)#img != 0
    queued = np.ascontiguousarray(np.where(img!=0, 0, -1)) #no voxel has been queued
    print(queued.shape)
    
    markedheap = []
    print("Now looping")
    start_time = time.time()
    innerloop_start_time = time.time()

    for (x,y,z) in zip(*nz): #iterate over nonzero and get redundant boundaries
        if isBoundary2(img, (x,y,z)):
            if isRedundant2(img, (x,y,z)):
                heapq.heappush(markedheap, (weightedImg[x,y,z],(x,y,z))) #add to heap with priority
                queued[x,y,z] = 1 #mark as queued

    print("--- selecting redundant boundary %s seconds ---" % (time.time() - innerloop_start_time))
    outerloop_start_time = time.time()
    while markedheap: #heap not empty
        t = heapq.heappop(markedheap)[-1]
        
        if isRedundant2(img, t):
            #check if is endpoint
            img[t] = 0 #deleting
            indexes26 = getNNeighIndexes2(t, 26, False)
            for p in indexes26:
                if queued[p] == 0: #is foreground but not queued
                    if isRedundant2(img, p):
                        heapq.heappush(markedheap, (weightedImg[p],p))
                        queued[p] = 1 #mark as queued
    print("--- the other pass %s seconds ---" % (time.time() - outerloop_start_time))
    #print(img)
    print("--- the complete loop took %s seconds ---" % (time.time() - start_time))


def newOrderedThinning4(img, weightedImg): #something broken here
    '''Ordered thinning algorithm considering the weightedImg values as order priority.
    Loop implemented according to the Larrabide 2007 paper.
    Optimised with heapq and some function calling reduction'''
    nz = np.nonzero(img)#img != 0
    queued = np.ascontiguousarray(np.where(img!=0, 0, -1)) #no voxel has been queued
    
    markedheap = []
    print("Now looping")
    start_time = time.time()
    innerloop_start_time = time.time()

    for (x,y,z) in zip(*nz): #iterate over nonzero and get redundant boundaries
        if isBoundary2(img, (x,y,z)):
            if isRedundant2(img, (x,y,z)):
                heapq.heappush(markedheap, (weightedImg[x,y,z],(x,y,z))) #add to heap with priority
                queued[x,y,z] = 1 #mark as queued

    print("--- selecting redundant boundary %s seconds ---" % (time.time() - innerloop_start_time))
    outerloop_start_time = time.time()
    while markedheap: #heap not empty
        t = heapq.heappop(markedheap)[-1]
        
        if isRedundant2(img, t): #we're checking again because maybe now some neighbour changed and is not redundant anymore
            #TO DO: check if is endpoint
            img[t] = 0 #deleting
            #check neighbours that might be redundant too
            for i in range(-1,2): 
                for j in range(-1,2):
                    for k in range(-1,2):
                        if (abs(i) + abs(j) + abs(k)) != 0: #hardcoding 26-neigh to avoid function calling?
                            if queued[t[0]+i,t[1]+j,t[2]+k] == 0: #is foreground but not queued
                                if isRedundant2(img, (t[0]+i,t[1]+j,t[2]+k)):
                                    heapq.heappush(markedheap, (weightedImg[t[0]+i,t[1]+j,t[2]+k],(t[0]+i,t[1]+j,t[2]+k)))
                                    queued[t[0]+i,t[1]+j,t[2]+k] = 1 #mark as queued
    print("--- the other pass %s seconds ---" % (time.time() - outerloop_start_time))
    #print(img)
    print("--- the complete loop took %s seconds ---" % (time.time() - start_time))


def newOrderedThinning5(img, weightedImg):
    '''Ordered thinning algorithm considering the weightedImg values as order priority.
    Loop implemented according to the Larrabide 2007 paper.
    Optimised with heapq and some function calling reduction'''
    nz = np.nonzero(img)#img != 0
    queued = np.ascontiguousarray(np.where(img!=0, 0, -1)) #no voxel has been queued
    print(queued.shape)
    
    markedheap = []
    print("Now looping")
    start_time = time.time()
    innerloop_start_time = time.time()

    for (x,y,z) in zip(*nz): #iterate over nonzero and get redundant boundaries
        if isBoundary2(img, (x,y,z)):
            if isRedundant2(img, (x,y,z)):
                heapq.heappush(markedheap, (weightedImg[x,y,z],(x,y,z))) #add to heap with priority
                queued[x,y,z] = 1 #mark as queued

    print("--- selecting redundant boundary %s seconds ---" % (time.time() - innerloop_start_time))
    outerloop_start_time = time.time()
    while markedheap: #heap not empty
        t = heapq.heappop(markedheap)[-1]
        
        if isRedundant2(img, t):
            #check if is endpoint
            img[t] = 0 #deleting
            for i in range(-1,2): 
                for j in range(-1,2):
                    for k in range(-1,2):
                        if (abs(i) + abs(j) + abs(k)) != 0: #hardcoding 26-neigh to avoid function calling?
                            if queued[t[0]+i,t[1]+j,t[2]+k] == 0: #is foreground but not queued
                                if isRedundant2(img, (t[0]+i,t[1]+j,t[2]+k)):
                                    heapq.heappush(markedheap, (weightedImg[t[0]+i,t[1]+j,t[2]+k],(t[0]+i,t[1]+j,t[2]+k)))
                                    queued[t[0]+i,t[1]+j,t[2]+k] = 1 #mark as queued
            
    print("--- the other pass %s seconds ---" % (time.time() - outerloop_start_time))
    #print(img)
    print("--- the complete loop took %s seconds ---" % (time.time() - start_time))


def newOrderedThinning6(img, weightedImg):


    def isBoundary2(img, v):
        '''Check if the 26-neighbourhood of v has at least on background voxel'''
        for i in range(-1, 2): 
            for j in range(-1, 2):
                for k in range(-1, 2):
                    if img[v[0]+i,v[1]+j,v[2]+k] == 0:
                        return True
        return False


    def isBackgroundSimple2(img, v):
        '''Check if the voxels of *value* in the *neigh* neighbourhood of *v* voxel are a connected component'''

        indexes18 = []
        indexes6 = 0
        #hardcoding neighbourhood extraction to avoid function calling
        for i in range(-1,2): 
            for j in range(-1,2):
                for k in range(-1,2):
                    if (abs(i) + abs(j) + abs(k)) != 0: #excluding center v
                        if (abs(i) + abs(j) + abs(k)) <=2: #can only move in two directions at a time for 18-neigh
                            if img[v[0]+i,v[1]+j,v[2]+k] == 0:
                                indexes18.append((v[0]+i,v[1]+j,v[2]+k)) #we only care about 18-neigh that's background
                                if (abs(i) + abs(j) + abs(k)) <= 1: #can only move in one direction at a time for 6-neigh
                                    indexes6 = indexes6 + 1
        
        if indexes6 == 0: #there are no background in the 6 connected
            return False
        
        else: #checked that part, now check if 18-neigh background is connected (flood fill)
            visit = [indexes18.pop()]
            while len(visit) != 0: #go over all the connected neighbours
                p = visit.pop()
                #print(p)
                for e in indexes18[:]: #check if any of the not visited is a neighbour (iterate over a copy so we can delete elements)
                    #print(e)
                    if max(abs(p[0]-e[0]), abs(p[1]-e[1]), abs(p[2]-e[2])) == 1:
                        visit.append(e)
                        indexes18.remove(e)
            if len(indexes18) == 0: #I visited all the background voxels in one connected pass
                return True
            else:
                return False


    def isBackgroundSimple3(img, v):
        '''Check if the voxels of *value* in the *neigh* neighbourhood of *v* voxel are a connected component'''

        indexes18 = []
        indexes6 = 0
        #hardcoding neighbourhood extraction to avoid function calling
        for i in range(-1,2): 
            for j in range(-1,2):
                for k in range(-1,2):
                    if (abs(i) + abs(j) + abs(k)) != 0: #excluding center v
                        if (abs(i) + abs(j) + abs(k)) <=2: #can only move in two directions at a time for 18-neigh
                            if img[v[0]+i,v[1]+j,v[2]+k] == 0:
                                indexes18.append((v[0]+i,v[1]+j,v[2]+k)) #we only care about 18-neigh that's background
                                if (abs(i) + abs(j) + abs(k)) <= 1: #can only move in one direction at a time for 6-neigh
                                    indexes6 = indexes6 + 1
        
        if indexes6 == 0: #there are no background in the 6 connected
            return False
        
        else: #checked that part, now check if 18-neigh background is connected (flood fill)
            
            visit = []
            heapq.heappush(visit, indexes18.pop()) #add to heap with priority[]
            while visit: #go over all the connected neighbours
                p = heapq.heappop(visit)
                delete = []
                for e in indexes18: #check if any of the not visited is a neighbour (iterate over a copy so we can delete elements)
                    #print(e)
                    if (abs(p[0]-e[0]) + abs(p[1]-e[1]) + abs(p[2]-e[2])) == 1: #moving on 6-connected manner
                        heapq.heappush(visit, e)
                        heapq.heappush(delete, e)
                while delete:
                    indexes18.remove(heapq.heappop(delete))
            if len(indexes18) == 0: #I visited all the background voxels in one connected pass
                return True
            else:
                return False



    def isForegroundSimple2(img, v):
        '''Check if the voxels of foreground in the 26 neighbourhood of *v* voxel are a connected component
        (which means, the removal of v does not affect the foreground connectivity)'''
        numneigh=0
        indexes=[]
        for i in range(-1, 2): 
            for j in range(-1, 2):
                for k in range(-1, 2):
                    if img[v[0]+i,v[1]+j,v[2]+k] == 1:
                        numneigh=numneigh+1
                        if (abs(i) + abs(j) + abs(k)) != 0: #we're excluding the center point v
                            indexes.append((v[0]+i,v[1]+j,v[2]+k))

        #if no foreground, v is not simple. If only one additional foreground, v is endline. If all neigh are foreground, v is not simple.
        if numneigh<=2 or numneigh == 27: 
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


    def isForegroundSimple3(img, v):
        '''Check if the voxels of foreground in the 26 neighbourhood of *v* voxel are a connected component
        (which means, the removal of v does not affect the foreground connectivity)'''
        numneigh=0
        indexes=[]
        for i in range(-1, 2): 
            for j in range(-1, 2):
                for k in range(-1, 2):
                    if img[v[0]+i,v[1]+j,v[2]+k] == 1:
                        numneigh=numneigh+1
                        if (abs(i) + abs(j) + abs(k)) != 0: #we're excluding the center point v
                            indexes.append((v[0]+i,v[1]+j,v[2]+k))

        #if no foreground, v is not simple. If only one additional foreground, v is endline. If all neigh are foreground, v is not simple.
        if numneigh<=2 or numneigh == 27: 
            return False
        
        visit = []
        heapq.heappush(visit, indexes.pop()) #pick one voxel to study connectivity
        while visit: #go over all the connected neighbours (flood fill)
            p = heapq.heappop(visit)
            #print(p)
            for e in indexes: #check if any of the not visited is a neighbour (iterate over a copy so we can delete elements in indexes)
                #print(e)
                delete = []
                if max(abs(p[0]-e[0]),abs(p[1]-e[1]),abs(p[2]-e[2])) == 1:
                    heapq.heappush(visit, e)
                    heapq.heappush(delete, e)
                while delete:
                    indexes.remove(heapq.heappop(delete))
        if len(indexes) == 0: #I visited all the foreground voxels in one connected pass
            return True
        else:
            return False




    def isEndPoint2(img, v):
        '''Check if the voxels of foreground in the 26 neighbourhood of *v* voxel are a connected component
        (which means, the removal of v does not affect the foreground connectivity)'''
        numneigh=0
        for i in range(-1, 2): 
            for j in range(-1, 2):
                for k in range(-1, 2):
                    if img[v[0]+i,v[1]+j,v[2]+k] != 0:
                        numneigh = numneigh+1

        #If only one additional foreground, v is endline.
        if numneigh<=2: 
            return True
        return False


    def isRedundant2(img, v):
        if isForegroundSimple2(img, v) and isBackgroundSimple2(img, v):
            return True
        else:
            return False

    '''Ordered thinning algorithm considering the weightedImg values as order priority.
    Loop implemented according to the Larrabide 2007 paper.
    Optimised with heapq and some function calling reduction'''
    nz = np.nonzero(img)#img != 0
    queued = np.ascontiguousarray(np.where(img!=0, 0, -1)) #no voxel has been queued
    print(queued.shape)
    
    markedheap = []
    print("Now looping")
    start_time = time.time()
    innerloop_start_time = time.time()

    for (x,y,z) in zip(*nz): #iterate over nonzero and get redundant boundaries
        if isBoundary2(img, (x,y,z)):
            if isForegroundSimple2(img, (x,y,z)) and isBackgroundSimple2(img, (x,y,z)):
                heapq.heappush(markedheap, (weightedImg[x,y,z],(x,y,z))) #add to heap with priority
                queued[x,y,z] = 1 #mark as queued

    #print("--- selecting redundant boundary %s seconds ---" % (time.time() - innerloop_start_time))
    outerloop_start_time = time.time()
    while markedheap: #heap not empty
        t = heapq.heappop(markedheap)[-1]
        
        if isForegroundSimple2(img, t) and isBackgroundSimple2(img, t):
            #if not isEndPoint2(img, t):
            img[t] = 0 #deleting
            for i in range(-1,2): 
                for j in range(-1,2):
                    for k in range(-1,2):
                        if (abs(i) + abs(j) + abs(k)) != 0: #hardcoding 26-neigh to avoid function calling?
                            if queued[t[0]+i,t[1]+j,t[2]+k] == 0: #is foreground but not queued
                                if isForegroundSimple2(img, (t[0]+i,t[1]+j,t[2]+k)) and isBackgroundSimple2(img, (t[0]+i,t[1]+j,t[2]+k)):
                                    heapq.heappush(markedheap, (weightedImg[t[0]+i,t[1]+j,t[2]+k],(t[0]+i,t[1]+j,t[2]+k)))
                                    queued[t[0]+i,t[1]+j,t[2]+k] = 1 #mark as queued
            
    #print("--- the other pass %s seconds ---" % (time.time() - outerloop_start_time))
    #print(img)
    print("--- the complete loop took %s seconds ---" % (time.time() - start_time))



def getWeightedImageEuclideanDistanceTransform(img):
    '''Generate a new image where the foreground voxels are weighted according to the squared euclidean distance
    (relevant for their ordered consideration to perform a binary thinning)'''
    profiled = edt.edtsq(np.ascontiguousarray(img), black_border=True, parallel=1, order="K")
    return profiled


def getWeightedImageEuclidean(img, profiled):
    '''Generate a new image where the foreground voxels are weighted according to the squared euclidean distance
    (relevant for their ordered consideration to perform a binary thinning)'''
    specified = np.nonzero(img) #get nonzero voxels
    for (x,y,z) in zip(*specified): #iterate over nonzero voxels
        profiled[x,y,z] = getMinEuclideanToBackground(img, (x,y,z))
        #profiled[x,y,z] = getProfileMeasure(img, (x,y,z), getMinEuclideanToBackground(img, (x,y,z)))

def getWeightedImageGrayscale(img, original):
    '''Generate a new image where the foreground voxels are assigned their (original) grayscale value
    (relevant for their ordered consideration to perform a binary thinning)'''
    weighted = np.zeros(img.shape, dtype=np.intc)
    specified = np.nonzero(img) #get nonzero voxels in segmentation
    for (x,y,z) in zip(*specified): #iterate over nonzero voxels
        weighted[x,y,z] = original[x,y,z]
    return weighted

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
    parser.add_argument("-rad", help="repeat closing", default=2, type=int, required=False)
    parser.add_argument("--gaussian", help="perform gaussian smoothing", action="store_true", required=False, default=False)
    #parser.add_argument("--connected", help="extract connected components", action="store_true", required=False, default=False)


    args = parser.parse_args()
    inputf = args.ifile
    outputf = args.ofile
    closing = args.closing
    inputgrayscale = args.gsfile
    rp = args.rp
    rad = args.rad
    gaussian = args.gaussian
    #conncomp = args.connected
    print("Starting thinning")
    img = readNIFTIasSITK(inputf)
    npimgorig = SITKToNumpy(readNIFTIasSITK(inputgrayscale)) if inputgrayscale is not None else None
    if closing:
        for _ in range(rp):
            img = binClosing(img, int(rad))
    if gaussian:
        img = gaussianSmoothDiscrete(img)
        npimg = SITKToNumpy(img)
        npthinned = binarySegmToBinarySkeleton3(npimg, npimgorig, 0.25)
    else:
        npimg = SITKToNumpy(img)
        npthinned = binarySegmToBinarySkeleton3(npimg, npimgorig)

    thinnedsitk = numpyToSITK(npthinned)
    thinnedsitk.SetOrigin(img.GetOrigin())
    thinnedsitk.SetSpacing(img.GetSpacing())
    thinnedsitk.SetDirection(img.GetDirection())

    writeSITK(thinnedsitk,outputf)


if __name__ == "__main__":
    main()

