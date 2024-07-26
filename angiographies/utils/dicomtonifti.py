import SimpleITK as sitk
import argparse
import os
from collections import OrderedDict
import json
import hashlib
from datetime import datetime

import os
import re

def hasDigits(inputString):
    return bool(re.search(r'\d', inputString))

join = os.path.join

def makedir(outputpath):
    os.makedirs(outputpath, exist_ok = True) #if the folder doesn't exist, we make it

def isdir(dirpath):
    return os.path.isdir(dirpath)

def isfile(filepath):
    return os.path.isfile(filepath)

def getsubdirs(ipath, join = True, contains = None, prefix = None, suffix = None, number = False, sort = True):
    '''join: return dir merged with path
    contains: check that this is part of the path
    prefix: check if this is at the beginning of the path
    suffix: check if this is at the end of the path
    number: check if dir has digit
    sort: return dirs sorted'''
    
    if join:
        l = os.path.join
    else:
        l = lambda x, y: y
    res = [l(ipath, i) for i in os.listdir(ipath) if os.path.isdir(os.path.join(ipath, i))
            and (prefix is None or i.startswith(prefix))
            and (suffix is None or i.endswith(suffix))
            and (contains is None or contains in i)
            and (not number or hasDigits(i))]
    if sort:
        res.sort()
    return res

def getsubfiles(ipath, join=True, contains = None, prefix=None, suffix=None, sort=True):
    '''join: return filename merged with path
    contains: check that this is part of the path
    prefix: check if this is at the beginning of the path
    suffix: check if this is at the end of the path
    sort: return files sorted'''

    if join:
        l = os.path.join
    else:
        l = lambda x, y: y
    res = [l(ipath, i) for i in os.listdir(ipath) if os.path.isfile(os.path.join(ipath, i))
            and (prefix is None or i.startswith(prefix))
            and (suffix is None or i.endswith(suffix))
            and (contains is None or contains in i)]
    if sort:
        res.sort()
    return res


def main():
    '''This script receives the path to a folder and recursively finds all dicom files to save as niftii'''

    parser = argparse.ArgumentParser()
    parser.add_argument("-idir", help="path to root", default="", required=True)
    #parser.add_argument("-ifilegraph", help="path to graph to convert to skeleton", default=None, required=False)
    parser.add_argument("-odir", help="name to output json", default=None, required=False)
    #parser.add_argument("-special", help="which folder are we segmenting", default=None, required=False)
    #parser.add_argument("-contain", help="seek for something in path", default=None, required=False)


    args = parser.parse_args()
    idir = args.idir
    odir = args.odir

    for root, subd, files in os.walk(idir): #recursively walk over the directory
        if len(subd)==0 and len(files)!=0:# and files[0][-3:]=="dcm" and not "2D" in root:# and "PRE" in root:

            print("read: ", root)
            reader = sitk.ImageSeriesReader()
            dicom_names = reader.GetGDCMSeriesFileNames(root)
            reader.SetFileNames(dicom_names)
            reader.MetaDataDictionaryArrayUpdateOn()
            #reader.LoadPrivateTagsOn()
            reader.Execute()

            # Get the sorted file names, opens all files in the directory and reads the meta-information
            # without reading the bulk pixel data
            series_ID = reader.GetMetaData(0,'0020|000e')
            sorted_file_names = sitk.ImageSeriesReader.GetGDCMSeriesFileNames(root, series_ID)
            print(len(sorted_file_names))

            # Read the bulk pixel data
            img = sitk.ReadImage(sorted_file_names)

            #reader = sitk.ImageSeriesReader()
            #dicom_names = reader.GetGDCMSeriesFileNames(root)
            #reader.SetFileNames(dicom_names)
            #reader.Execute()
            # Get the sorted file names
            #sorted_file_names = sitk.ImageSeriesReader.GetGDCMSeriesFileNames(root, series_ID)
            #series_reader = sitk.ImageSeriesReader()
            #series_reader.SetFileNames(sorted_file_names)
            #series_reader.MetaDataDictionaryArrayUpdateOn()
            #series_reader.LoadPrivateTagsOff()
            #image_dicom = reader.Execute()
            print("size: ", img.GetSize())
            print("spacing: ", img.GetSpacing())
            patient = root.split(idir)[-1].split(os.path.sep)[1].split(" ")[0]
            #treatment = root.split(idir)[-1].split(os.path.sep)[2]
            print(os.path.join(odir,patient+"-post.nii.gz")) 
            sitk.WriteImage(img, os.path.join(odir,patient+"-post.nii.gz"))
    
if __name__ == "__main__":
    main()