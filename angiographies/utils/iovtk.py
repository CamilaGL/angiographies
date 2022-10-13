import vtk
from vtkmodules.util import numpy_support
import numpy

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
    #print(reader.GetNIFTIHeader().GetQOffsetX())
    #print(reader.GetQFormMatrix())
    #print(reader.GetQFormMatrix().GetElement(0,3),reader.GetQFormMatrix().GetElement(1,3),reader.GetQFormMatrix().GetElement(2,3))
    inputData = reader.GetOutput()
    #remember vtk always returns 0 0 0 as origin
    return inputData


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

def writeNIFTIwithVTK(image, fileName):
    '''TO DO: test!'''
    #image as vtkimage
    writer = vtk.vtkNIFTIImageWriter()
    writer.SetInputData(image)
    writer.SetFileName(fileName)
    writer.Write()

