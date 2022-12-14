'''
This script segments a volume with a simple threshold.

requires:
    vtk

'''
import vtk
import os
import argparse
import os

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-ipath", help="path to file", default="", required=True)
    parser.add_argument("-opath", help="path to output file", default="", required=True)
    parser.add_argument("-thresh", help="threshold value (keep equal or greater)", type=int, default=1000, required=False)

    args = parser.parse_args()
    ipath = args.ipath
    opath = args.opath
    thresh = args.thresh

    if os.path.isfile(ipath):
        
        #read niftii
        reader = vtk.vtkNIFTIImageReader()
        reader.SetFileName(ipath)
        reader.Update()
        print(reader.GetOutput().GetDimensions())

        #threshold
        threshold = vtk.vtkImageThreshold()
        threshold.SetInputConnection(reader.GetOutputPort())
        threshold.ThresholdByUpper(thresh)  #th
        threshold.ReplaceInOn()
        threshold.SetInValue(1)  # set all values above or equal th to 1
        threshold.ReplaceOutOn()
        threshold.SetOutValue(0)  # set all values below th to 1
        threshold.Update()

        #write as nifti
        writer = vtk.vtkNIFTIImageWriter()
        writer.SetInputConnection(threshold.GetOutputPort())
        writer.SetFileName(opath)
        writer.Write()


if __name__ == "__main__":
    main()
