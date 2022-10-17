import SimpleITK as sitk
import numpy
import math
from angiographies.utils.formatconversionsitk import SITKToNumpy



def binClosing(binImg, rad):
    vectorRadius=(rad,rad,rad)
    kernel=sitk.sitkBall
    print(binImg.GetPixelIDTypeAsString())
    if ("integer" not in binImg.GetPixelIDTypeAsString()):
        return sitk.BinaryMorphologicalClosing(sitk.Cast(binImg, sitk.sitkInt32),vectorRadius, kernel)
    return sitk.BinaryMorphologicalClosing(binImg,vectorRadius, kernel)

def thresholdImageSITK(inputImage, lt, ut):
    #input as itk
    tf = sitk.ThresholdImageFilter()
    #tf.SetOutsideValue(-4000)
    tf.SetLower(lt)
    tf.SetUpper(ut)
    segm1 = tf.Execute(inputImage) #all background is black
    tf2 = sitk.ThresholdImageFilter()
    tf.SetOutsideValue(1) #not background is 1
    tf.SetLower(0)
    tf.SetUpper(0)
    #tf.ThresholdAbove()
    segm = tf.Execute(segm1) #all values not background are 1
    #tf.SetOutsideValue(0) #background to 0
    #tf.SetLower(1)
    #tf.SetUpper(1) #everything not foreground is background
    #segm2 = tf.Execute(segm)
    #vectorRadius=(1,1,1)
    #kernel=sitk.sitkBall
    # seg_implicit_thresholds_clean = sitk.BinaryMorphologicalClosing(segm, 
    #                                                                 vectorRadius,
    #                                                                 kernel)
    #print(segm.dtype)
    return segm #seg_implicit_thresholds_clean

def getLargestConnected(inputImage):
    #inputimage as sitk (change to numpy?)
    #is this efficient? absolutely not, but at least it works.
    connectedFilter = sitk.ConnectedComponentImageFilter()
    inputImageconn = connectedFilter.Execute(inputImage)
    print(connectedFilter.GetObjectCount())
    connnp = SITKToNumpy(inputImageconn)
    val, count = numpy.unique(connnp, return_counts=True)
    values = zip(val, count)
    v = 0
    c = 0
    for i, j in values:
        if j>c and i!=0:
            c=j
            v=int(i)
    print(v)
    return thresholdImageSITK(inputImageconn, v, v)


def gaussianSmoothRecursive(binImg, m=1.7):
    spa = binImg.GetSpacing()[0]
    gaussian = sitk.SmoothingRecursiveGaussianImageFilter()
    gaussian.SetSigma(spa*m)
    blurredImage = gaussian.Execute(binImg)
    return blurredImage


def gaussianSmoothDiscrete(binImg, m=1.7):
    #spa = binImg.GetSpacing()[0]
    gaussian = sitk.DiscreteGaussianImageFilter()
    gaussian.SetVariance(math.sqrt(m))
    gaussian.SetMaximumKernelWidth(1)
    gaussian.SetUseImageSpacing(True)
    blurredImage = gaussian.Execute(sitk.Cast(binImg, sitk.sitkFloat32))
    return blurredImage
