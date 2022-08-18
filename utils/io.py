import vtk
import SimpleITK as sitk
from vtk.util import numpy_support
import numpy
import json

def readVTKasVTK(fileName):
    #read with VTK file as VTK Image
    dreader = vtk.vtkGenericDataObjectReader()
    dreader.SetFileName(fileName)
    dreader.Update()

    imageVTK= dreader.GetOutput()
    dimensions = imageVTK.GetDimensions()
    spacing = imageVTK.GetSpacing()
    origin = imageVTK.GetOrigin()

    print(dimensions)
    print(spacing)
    print(origin)
    return imageVTK

def readSTLasVTK(filename):
    dreader = vtk.vtkSTLReader()
    dreader.SetFileName(filename)
    dreader.Update()
    return dreader.GetOutput()

def readDICOMasSITK(inputDirectory):
    #TO DO
    reader = sitk.ImageSeriesReader()
    dicom_names = reader.GetGDCMSeriesFileNames(inputDirectory)
    reader.SetFileNames(dicom_names)
    reader.ReadImageInformation()

    # Get the sorted file names, opens all files in the directory and reads the meta-information
    # without reading the bulk pixel data
    series_ID = reader.GetMetaData('0020|000e')
    sorted_file_names = sitk.ImageSeriesReader.GetGDCMSeriesFileNames(inputDirectory, series_ID)

    # Read the bulk pixel data
    img = sitk.ReadImage(sorted_file_names)
    # reader.GetMetaData(0, '0010|0010')

def readVTIasVTKFlipped(inputDirectory):
    #read DICOM as VTK Image, correcting z order
    dreader = vtk.vtkXMLImageDataReader()
    dreader.SetFileName(inputDirectory)
    dreader.UpdateWholeExtent()
    dreader.Update()

    imageVTK= dreader.GetOutput()
    dimensions = imageVTK.GetDimensions()
    spacing = imageVTK.GetSpacing()
    origin = imageVTK.GetOrigin()

    print(dimensions)
    print(spacing)
    print(origin)

    data = imageVTK.GetPointData() # get VTK-formatted data
    data = numpy_support.vtk_to_numpy(data.GetArray(0)) # convert VTKdata to numpy data
    print(data.shape)
    
    dataReshaped = data.reshape((dimensions[2],dimensions[1],dimensions[0])) #reshape data from 1D -> 3D 
    dataReordered = numpy.ascontiguousarray(dataReshaped[:, :, :])
    dataRavelled = dataReordered.ravel()
    flippedDataArray = numpy_support.numpy_to_vtk(dataRavelled, True)
    flippedData = vtk.vtkImageData()
    flippedData.GetPointData().SetScalars(flippedDataArray)
    flippedData.SetSpacing(spacing)
    flippedData.SetOrigin(origin)
    flippedData.SetDimensions(dimensions)

    return flippedData


def readDICOMasVTKFlipped(inputDirectory):
    #read DICOM as VTK Image, correcting z order
    dreader = vtk.vtkDICOMImageReader()
    dreader.SetDirectoryName(inputDirectory)
    dreader.Update()

    imageVTK= dreader.GetOutput()
    dimensions = imageVTK.GetDimensions()
    spacing = imageVTK.GetSpacing()
    origin = imageVTK.GetOrigin()

    print(dimensions)
    print(spacing)
    print(origin)

    data = imageVTK.GetPointData() # get VTK-formatted data
    data = numpy_support.vtk_to_numpy(data.GetArray(0)) # convert VTKdata to numpy data
    print(data.shape)
    
    dataReshaped = data.reshape((dimensions[2],dimensions[1],dimensions[0])) #reshape data from 1D -> 3D 
    dataReordered = numpy.ascontiguousarray(dataReshaped[::-1, :, :])
    dataRavelled = dataReordered.ravel()
    flippedDataArray = numpy_support.numpy_to_vtk(dataRavelled, True)
    flippedData = vtk.vtkImageData()
    flippedData.GetPointData().SetScalars(flippedDataArray)
    flippedData.SetSpacing(spacing)
    flippedData.SetOrigin(origin)
    flippedData.SetDimensions(dimensions)

    return flippedData


def readVTPPolydata(fileName):
    #read with VTK file as VTK Image
    dreader = vtk.vtkXMLPolyDataReader()
    dreader.SetFileName(fileName)
    dreader.Update()
    #dwriter.WriteExtentOn()
    #dwriter.SetWholeExtent(inputImage.GetWholeExtent())
    inputData = dreader.GetOutput()
    return inputData


def readNIFTIasVTK(fileName):
    reader = vtk.vtkNIFTIImageReader()
    reader.SetFileName(fileName)
    reader.Update()
    inputData = reader.GetOutput()
    return inputData


def readNIFTIasSITK(filename):
    return sitk.ReadImage(filename)


def writeVTKPolydataasVTK(inputImage, fileName):
    #read with VTK file as VTK Image
    dwriter = vtk.vtkPolyDataWriter()
    dwriter.SetFileName(fileName)
    #dwriter.WriteExtentOn()
    #dwriter.SetWholeExtent(inputImage.GetWholeExtent())
    dwriter.SetInputData(inputImage)
    dwriter.Write()


def writeVTKPolydataasVTP(inputImage, fileName):
    #read with VTK file as VTK Image
    dwriter = vtk.vtkXMLPolyDataWriter()
    dwriter.SetFileName(fileName)
    #dwriter.WriteExtentOn()
    #dwriter.SetWholeExtent(inputImage.GetWholeExtent())
    dwriter.SetInputData(inputImage)
    dwriter.Write()

def writeSITK(imageSITK, outputFile):
    sitk.WriteImage(imageSITK, outputFile)

def writejson(dict, outputFile):
    with open(outputFile, 'w') as f:
        json.dump(dict, f, sort_keys=True, indent=4)

def readjson(inputFile):
    with open(inputFile, 'r') as f:
        a = json.load(f)
    return a