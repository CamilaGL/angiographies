'''
This script receives either a folder or a file path, to perform screen captures of 3D volumes and save them to the output folder path.

requires:
    numpy
    vtk
    PIL
    moviepy

'''
import vtk
import os
import argparse
from fileoperations import *

def tupleToString(tup):
    st = ''.join(map(str, tup))
    return st

# def screencap(ifile, opath):
#     '''from https://gist.github.com/pangyuteng/facd430d0d9761fc67fff4ff2e5fffc3 '''

#     reader = vtk.vtkNIFTIImageReader()
#     reader.SetFileName(ifile)
#     reader.Update()

#     threshold = vtk.vtkImageThreshold()
#     threshold.SetInputConnection(reader.GetOutputPort())
#     threshold.ThresholdByLower(50)  #th
#     threshold.ReplaceInOn()
#     threshold.SetInValue(0)  # set all values below th to 0
#     threshold.ReplaceOutOn()
#     threshold.SetOutValue(1)  # set all values above th to 1
#     threshold.Update()

#     '''
#     voi = vtk.vtkExtractVOI()
#     voi.SetInputConnection(threshold.GetOutputPort()) 
#     voi.SetVOI(0,95, 0,95, 0,59)
#     voi.SetSampleRate(1,1,1)
#     #voi.SetSampleRate(3,3,3)
#     voi.Update()#necessary for GetScalarRange()
#     srange= voi.GetOutput().GetScalarRange()#needs Update() before!
#     print("Range", srange)
#     '''

#     dmc = vtk.vtkDiscreteMarchingCubes()
#     dmc.SetInputConnection(threshold.GetOutputPort())
#     #dmc.SetInputConnection(voi.GetOutputPort())
#     dmc.GenerateValues(1, 1, 1)
#     dmc.Update()

#     mapper = vtk.vtkPolyDataMapper()
#     mapper.SetInputConnection(dmc.GetOutputPort())
#     mapper.Update()

#     actor = vtk.vtkActor()
#     actor.SetMapper(mapper)
#     mapper.ScalarVisibilityOff()

#     actor.GetProperty().SetColor(salmon)
#     actor.GetProperty().SetOpacity(0.5)
#     actor.RotateX(90)

#     renderer = vtk.vtkRenderer()
#     renderWindow = vtk.vtkRenderWindow()
#     renderWindow.SetOffScreenRendering(1)
#     renderWindow.AddRenderer(renderer)

#     renderer.AddActor(actor)
#     renderer.SetBackground(1.0, 1.0, 1.0)

#     center = actor.GetCenter()

#     camera = renderer.MakeCamera()
#     camera.SetPosition(0,0,800)
#     camera.SetFocalPoint(center)

#     renderer.SetActiveCamera(camera)
#     renderWindow.Render()

#     focal_point = camera.GetFocalPoint()
#     view_up = camera.GetViewUp()
#     position = camera.GetPosition() 

#     axis = [0,0,0]
#     axis[0] = -1*camera.GetViewTransformMatrix().GetElement(0,0)
#     axis[1] = -1*camera.GetViewTransformMatrix().GetElement(0,1)
#     axis[2] = -1*camera.GetViewTransformMatrix().GetElement(0,2)

#     print(position,focal_point,view_up,)

#     print(camera.GetViewTransformMatrix())
#     print(camera.GetViewTransformMatrix().GetElement(0,0))
#     print(camera.GetViewTransformMatrix().GetElement(0,1))
#     print(camera.GetViewTransformMatrix().GetElement(0,2))

#     frame_list = []
#     for n,q in enumerate([10]*35):

#         transform = vtk.vtkTransform()
#         transform.Identity()

#         transform.Translate(*center)
#         transform.RotateWXYZ(q,view_up)
#         transform.RotateWXYZ(0,axis)
#         transform.Translate(*[-1*x for x in center])

#         new_position = [0,0,0]
#         new_focal_point = [0,0,0]
#         transform.TransformPoint(position,new_position)
#         transform.TransformPoint(focal_point,new_focal_point)

#         camera.SetPosition(new_position)
#         camera.SetFocalPoint(new_focal_point)

#         focal_point = camera.GetFocalPoint()
#         view_up = camera.GetViewUp()
#         position = camera.GetPosition() 

#         camera.OrthogonalizeViewUp();
#         renderer.ResetCameraClippingRange();
        
#         renderWindow.Render()
#         windowToImageFilter = vtk.vtkWindowToImageFilter()
#         windowToImageFilter.SetInputBufferTypeToRGBA()
#         windowToImageFilter.SetInput(renderWindow)
#         windowToImageFilter.Update()

#         writer = vtk.vtkPNGWriter()
        
#         fpath = "zisosurface{}.png".format(n)
#         writer.SetFileName(os.path.join(opath,fpath))
#         writer.SetInputConnection(windowToImageFilter.GetOutputPort())
#         writer.Write()

#         frame_list.append(fpath)
        
#     duration = 2
#     fps = 17.5
#     time_list = list(np.arange(0,duration,1./fps))
#     img_dict = {a:f for a,f in zip(time_list,frame_list)}

#     import PIL
#     from moviepy import editor
#     def make_frame(t):    
#         fpath= img_dict[t]
#         im = PIL.Image.open(os.path.join(opath,fpath))
#         ar = np.asarray(im)
#         return ar
        
#     gif_path = os.path.join(opath,'ani.gif')
#     clip = editor.VideoClip(make_frame, duration=duration)
#     clip.write_gif(gif_path, fps=fps)

        
#     '''
#     ref
#     https://www.programcreek.com/python/example/9480/vtk.vtkRenderWindows
#     https://pyscience.wordpress.com/2014/09/11/surface-extraction-creating-a-mesh-from-pixel-data-using-python-and-vtk/
#     https://stackoverflow.com/questions/43046798/how-to-rotate-a-vtk-scene-around-a-specific-point
#     https://gamedev.stackexchange.com/questions/20758/how-can-i-orbit-a-camera-about-its-target-point
#     '''


def screencap2(ifile, opath, filename="isosurf"):
    '''from https://gist.github.com/pangyuteng/facd430d0d9761fc67fff4ff2e5fffc3 '''

    reader = vtk.vtkNIFTIImageReader()
    reader.SetFileName(ifile)
    reader.Update()
    print(reader.GetOutput().GetDimensions())

    threshold = vtk.vtkImageThreshold()
    threshold.SetInputConnection(reader.GetOutputPort())
    threshold.ThresholdByUpper(1)  #th
    threshold.ReplaceInOn()
    threshold.SetInValue(1)  # set all values below th to 0
    threshold.ReplaceOutOn()
    threshold.SetOutValue(0)  # set all values above th to 1
    threshold.Update()

    '''
    voi = vtk.vtkExtractVOI()
    voi.SetInputConnection(threshold.GetOutputPort()) 
    voi.SetVOI(0,95, 0,95, 0,59)
    voi.SetSampleRate(1,1,1)
    #voi.SetSampleRate(3,3,3)
    voi.Update()#necessary for GetScalarRange()
    srange= voi.GetOutput().GetScalarRange()#needs Update() before!
    print("Range", srange)
    '''

    dmc = vtk.vtkDiscreteMarchingCubes()
    dmc.SetInputConnection(threshold.GetOutputPort())
    #dmc.SetInputConnection(voi.GetOutputPort())
    dmc.GenerateValues(1, 1, 1)
    dmc.Update()

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(dmc.GetOutputPort())
    mapper.Update()

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    mapper.ScalarVisibilityOff()

    actor.GetProperty().SetColor(1,0,0)
    actor.GetProperty().SetOpacity(0.8)
    actor.RotateX(90)

    renderer = vtk.vtkRenderer()
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.SetOffScreenRendering(1)
    renderWindow.AddRenderer(renderer)
    renderWindow.SetAlphaBitPlanes(1)
    renderWindow.SetSize(1500,1500)

    renderer.AddActor(actor)
    renderer.SetBackground(1.0, 1.0, 1.0)

    center = actor.GetCenter()

    camera = renderer.MakeCamera()
    #camera.SetPosition(0,0,800)
    camera.SetPosition(400,150,0)
    camera.SetFocalPoint(center)

    renderer.SetActiveCamera(camera)
    renderWindow.Render()
    focal_point = camera.GetFocalPoint()
    view_up = camera.GetViewUp()
    view_up = (0.0, -1.0, -1.0) 
    first_view_up = view_up
    camera.SetViewUp(view_up)
    position = camera.GetPosition() 

    axis = [0,0,0]
    axis[0] = -1*camera.GetViewTransformMatrix().GetElement(0,0)
    axis[1] = -1*camera.GetViewTransformMatrix().GetElement(0,1)
    axis[2] = -1*camera.GetViewTransformMatrix().GetElement(0,2)

    print(position,focal_point,view_up,)

    print(camera.GetViewTransformMatrix())
    print(camera.GetViewTransformMatrix().GetElement(0,0))
    print(camera.GetViewTransformMatrix().GetElement(0,1))
    print(camera.GetViewTransformMatrix().GetElement(0,2))

    frame_list = []
    for n,q in enumerate([60]*6):

        transform = vtk.vtkTransform()
        transform.Identity()

        transform.Translate(*center)
        transform.RotateWXYZ(q,view_up)
        transform.RotateWXYZ(0,axis)
        transform.Translate(*[-1*x for x in center])

        new_position = [0,0,0]
        new_focal_point = [0,0,0]
        transform.TransformPoint(position,new_position)
        transform.TransformPoint(focal_point,new_focal_point)

        camera.SetPosition(new_position)
        camera.SetFocalPoint(new_focal_point)

        focal_point = camera.GetFocalPoint()
        view_up = camera.GetViewUp()
        position = camera.GetPosition() 

        camera.OrthogonalizeViewUp();
        renderer.ResetCameraClippingRange();
        
        renderWindow.Render()
        windowToImageFilter = vtk.vtkWindowToImageFilter()
        #windowToImageFilter.SetInputBufferTypeToRGBA()
        windowToImageFilter.SetInput(renderWindow)
        windowToImageFilter.SetScale(1)
        windowToImageFilter.Update()

        writer = vtk.vtkPNGWriter()
        
        
        fpath = "{filename}_{frame}.png".format(frame = n, filename = filename)#, viewup = ''.join("{:+03.0f}".format(item*10.0) for item in first_view_up))
        writer.SetFileName(os.path.join(opath,fpath))
        writer.SetInputConnection(windowToImageFilter.GetOutputPort())
        writer.Write()

        frame_list.append(fpath)
        
    # duration = 2
    # fps = 17.5
    # time_list = list(np.arange(0,duration,1./fps))
    # img_dict = {a:f for a,f in zip(time_list,frame_list)}

    # import PIL
    # from moviepy import editor
    # def make_frame(t):    
    #     fpath= img_dict[t]
    #     im = PIL.Image.open(os.path.join(opath,fpath))
    #     ar = np.asarray(im)
    #     return ar
        
    # gif_path = os.path.join(opath,filename+".gif")
    # clip = editor.VideoClip(make_frame, duration=duration)
    # clip.write_gif(gif_path, fps=fps)

        
    '''
    ref
    https://www.programcreek.com/python/example/9480/vtk.vtkRenderWindows
    https://pyscience.wordpress.com/2014/09/11/surface-extraction-creating-a-mesh-from-pixel-data-using-python-and-vtk/
    https://stackoverflow.com/questions/43046798/how-to-rotate-a-vtk-scene-around-a-specific-point
    https://gamedev.stackexchange.com/questions/20758/how-can-i-orbit-a-camera-about-its-target-point
    '''



def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-ipath", help="path to file or folder to capture", default="", required=True)
    parser.add_argument("-opath", help="path to folder to save screencap", default="", required=True)
    parser.add_argument("--all", help="capture all files in folder", action="store_true", required=False, default=False)

    args = parser.parse_args()
    ipath = args.ipath
    opath = args.opath
    allfiles = args.all

    if allfiles and isdir(ipath):
        thefiles = getsubfiles(ipath, suffix=".nii.gz")
    elif isfile(ipath):
        thefiles = [ipath]
    else:
        print("Unable to read files")
        return

    for f in thefiles:
        filename = (f.split(os.path.sep)[-1]).split(".")[0]
        print(filename)
        print(f)
        print(opath)
        screencap2(f, opath, filename)



if __name__ == "__main__":
    main()
