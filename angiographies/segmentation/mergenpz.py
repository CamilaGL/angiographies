'''
This script merges npz files as an ensemble.

requires:
    numpy
    pickle
    sitk
'''

import numpy as np
import os
import argparse
from angiographies.utils.iopickle import *
from angiographies.utils.iositk import writeSITK
from angiographies.utils.formatconversionsitk import numpyToSITK
from angiographies.utils.fileoperations import *
import re

def hasDigits(inputString):
    return bool(re.search(r'\d', inputString))


def mergeSoftmax(files, picklefiles, ofile, overwrite, savenpz):
    if overwrite or not isfile(ofile):
        softmax = [np.load(f)['softmax'][None] for f in files]
        softmax = np.vstack(softmax)
        softmax = np.mean(softmax, 0) #mean across all cases
        pklfiles = [readPickle(f) for f in picklefiles]

        softmax = softmax.argmax(0) #because we have two dimensions representing two labels

        mergeditk = numpyToSITK(softmax.astype(np.uint8))
        mergeditk.SetSpacing(pklfiles[0]['itk_spacing'])
        mergeditk.SetOrigin(pklfiles[0]['itk_origin'])
        mergeditk.SetDirection(pklfiles[0]['itk_direction'])
        writeSITK(mergeditk, ofile)

        if savenpz:
            np.savez_compressed(ofile[:-7] + ".npz", softmax=softmax)
            writePickle(pklfiles, ofile[:-7] + ".pkl")


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-ifolder", help="path to folder with folds", default="", required=True)
    #parser.add_argument("-ofolder", help="path to folder with ensemble", default="", required=True)
    parser.add_argument("--overwrite", help="if exists, overwrite merged", action="store_true", required=False, default=False)

    args = parser.parse_args()
    ipath = args.ifolder
    overwrite = args.overwrite
    #outputf = args.ofolder

    print(ipath)
    if os.path.isdir(ipath): #check input path exists
        folds = getsubdirs(ipath, number=True) #get all folds
        print(folds)
        opath = join(ipath, "ensemble") #we're saving merged files here
        print(opath)
        makedir(opath) #if folder doesn't exist, makedir

        cases = [getsubfiles(i, suffix=".npz", join=False) for i in folds] #get cases in all folds
        cases = [i[:-4] for j in cases for i in j] #get all case numbers without file extension
        cases = np.unique(cases) #keep only one of each

        for i in cases: #for each case
            #check all folds have this case
            files = [join(f, i + ".npz") for f in folds if isfile(join(f, i + ".npz"))]
            picklefiles = [join(f, i + ".pkl") for f in folds if isfile(join(f, i + ".pkl"))]
            if (len(files) == len(picklefiles) == len(folds)): #all files exist in all folds
                print("este archivo", i, "esta para ser procesado")
                mergeSoftmax(files, picklefiles, join(opath,i+".nii.gz"), overwrite, savenpz=False)
            else:
                print("Not all npz or pkl are available in all folds for case "+i)



if __name__ == "__main__":
    main()
