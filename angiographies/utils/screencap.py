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

    actor.GetProperty().SetColor(0.74,0.08,0.08)
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
    #camera.SetPosition(400,150,0) # spacing 0.227
    camera.SetPosition(1600,600,0) # spacing 1
    camera.SetFocalPoint(center)

    renderer.SetActiveCamera(camera)
    renderWindow.Render()
    focal_point = camera.GetFocalPoint()
    view_up = camera.GetViewUp()
    view_up = (0.0, -1.0, -1.0) #0 1 1
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

def screencapMulticlassfg(ifile, output, casenum="137"):
    '''from https://gist.github.com/pangyuteng/facd430d0d9761fc67fff4ff2e5fffc3 
    Save screencap of imagenumpy in path to be used afterwards
    '''

    reader = vtk.vtkXMLImageDataReader()
    reader.SetFileName(ifile)
    reader.Update()

    colormapping = {"2": (0.10,0.3,0.65),
                "1": (0.9,0.8,0.7),
                "3": (1.0,0.4,0.15),
                "4": (0.11,0.82,0.68)}#(1.0,0.8,0.0)}
    opacitymapping = {"2": 0.8,
            "1": 0.8,
            "3": 0.8,
            "4": 0.8}


    renderer = vtk.vtkRenderer()
    renderWindow = vtk.vtkRenderWindow()
    rendererfg = vtk.vtkRenderer()
    renderWindow.SetNumberOfLayers(2)
    renderer.SetLayer(0)
    rendererfg.SetLayer(1)
    renderWindow.AddRenderer(renderer)
    renderWindow.AddRenderer(rendererfg)
    renderWindow.SetOffScreenRendering(1)
    renderWindow.SetAlphaBitPlanes(1)
    renderWindow.SetSize(2000,1500)
    imgrange = reader.GetOutput().GetScalarRange()

    actors = []
    center = (0,0,0)
    for u in range(int(imgrange[0])+1, int(imgrange[1])+1):#[1,2,3]:
        threshold = vtk.vtkImageThreshold()
        threshold.SetInputConnection(reader.GetOutputPort()) #Connection(reader.GetOutputPort())
        threshold.ThresholdBetween(u, u)  #th
        threshold.ReplaceInOn()
        threshold.SetInValue(1)  # set all values below th to 0
        threshold.ReplaceOutOn()
        threshold.SetOutValue(0)  # set all values above th to 1
        threshold.Update()

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

        actor.GetProperty().SetColor(colormapping[str(u)])

        actor.GetProperty().SetOpacity(opacitymapping[str(u)])

        #actor.RotateX(90)
        if u == 4:
            rendererfg.AddActor(actor)
        else:
            renderer.AddActor(actor)
    
    renderer.SetBackground(1.0, 1.0, 1.0)

    camera = renderer.MakeCamera()
    if casenum == "194":
        # camera.SetPosition(690.9347955520454, 1112.0667517889565, 626.2356060635063) # spacing 1
        # camera.SetFocalPoint(220.6152792316476, 312.97492897695184, 294.4309474200802) 
        # camera.SetViewUp(0.057999847002381484, -0.41176253845691324, 0.9094435824564487) 
        # camera.SetPosition(1016.3155743032796, 874.2863804281911, 351.6890125632976) # spacing 1
        # camera.SetFocalPoint(235.5091650595668, 279.1853777201952, 273.98014515282676) 
        # camera.SetViewUp(0.0042712830451648726, -0.13498874609218414, 0.9908379254800495)  
        # camera.SetPosition(176.89317125378793, 119.44684296667383, 338.50862258646816) # spacing 1
        # camera.SetFocalPoint(269.0349198391431, 253.34569168733125, 268.11962547264505) 
        # camera.SetViewUp(-0.2819577714595236, 0.5912390800225928, 0.7556031798289478)  
        # camera.SetPosition(-928.0086345059518, 120.61989189153255, 287.2433403515401) # spacing 1
        # camera.SetFocalPoint(255.26266497455592, 256.7569724646287, 251.43606504964842) 
        # camera.SetViewUp(0.04653638264257264, -0.1445083474366278, 0.9884086718618362)      
        camera.SetPosition(-0.07771454318191218, 211.83170928157543, 272.90595920292253) # spacing 1
        camera.SetFocalPoint(257.436595768047, 241.4591044756921, 265.11325282397195) 
        camera.SetViewUp(0.04653638264257264, -0.1445083474366278, 0.9884086718618362)                                         

    renderer.SetActiveCamera(camera)
    rendererfg.SetActiveCamera(camera)
    renderWindow.Render()

    camera.OrthogonalizeViewUp();
    rendererfg.ResetCameraClippingRange();
    renderer.ResetCameraClippingRange();
    
    renderWindow.Render()
    windowToImageFilter = vtk.vtkWindowToImageFilter()
    #windowToImageFilter.SetInputBufferTypeToRGBA()
    windowToImageFilter.SetInput(renderWindow)
    windowToImageFilter.SetScale(1)
    windowToImageFilter.Update()

    writer = vtk.vtkPNGWriter()
    writer.SetFileName(output)
    writer.SetInputConnection(windowToImageFilter.GetOutputPort())
    writer.Write()



def screencapMulticlass163(ifile, output, casenum="137"):
    '''from https://gist.github.com/pangyuteng/facd430d0d9761fc67fff4ff2e5fffc3 
    Save screencap of imagenumpy in path to be used afterwards
    '''

    reader = vtk.vtkXMLImageDataReader()
    reader.SetFileName(ifile)
    reader.Update()
    # print(reader.GetOutput().GetDimensions())
    #print(len(np.transpose(image.nonzero())))
    # colormapping = {"2": (0.74,0.08,0.08),
    #                 "1": (0.9,0.69,0.0),
    #                 "3": (0.0,0.45,0.5)}
    
        # colormapping = {"2": (0.15,0.32,0.45),
        #         "1": (0.9,0.8,0.7),
        #         "3": (1.0,0.4,0.15)}



    colormapping = {"2": (0.10,0.3,0.65),
                "1": (0.9,0.8,0.7),
                "3": (1.0,0.4,0.15)}
    opacitymapping = {"2": 0.8,
            "1": 0.8,
            "3": 1.0}


    renderer = vtk.vtkRenderer()
    renderWindow = vtk.vtkRenderWindow()
    rendererfg = vtk.vtkRenderer()
    renderWindow.SetNumberOfLayers(2)
    renderer.SetLayer(0)
    rendererfg.SetLayer(1)
    renderWindow.AddRenderer(renderer)
    renderWindow.AddRenderer(rendererfg)
    renderWindow.SetOffScreenRendering(1)
    renderWindow.SetAlphaBitPlanes(1)
    renderWindow.SetSize(2000,1500)

    actors = []
    center = (0,0,0)
    for u in [1,2,3]:
        threshold = vtk.vtkImageThreshold()
        threshold.SetInputConnection(reader.GetOutputPort()) #Connection(reader.GetOutputPort())
        threshold.ThresholdBetween(u, u)  #th
        threshold.ReplaceInOn()
        threshold.SetInValue(1)  # set all values below th to 0
        threshold.ReplaceOutOn()
        threshold.SetOutValue(0)  # set all values above th to 1
        threshold.Update()

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

        actor.GetProperty().SetColor(colormapping[str(u)])

        actor.GetProperty().SetOpacity(0.8)

        #actor.RotateX(90)
        if u == 2:
            rendererfg.AddActor(actor)
        else:
            renderer.AddActor(actor)
    
    renderer.SetBackground(1.0, 1.0, 1.0)

    camera = renderer.MakeCamera()

    camera.SetPosition(1082.843131000895, -152.90436517958975, -459.53254983613095) # spacing 1
    camera.SetFocalPoint(272.120231576054, 248.165519921031, 276.7944973668337) 
    camera.SetViewUp(0.2856371298967383, -0.6737763549117353, 0.6814960407707511)                                       

    renderer.SetActiveCamera(camera)
    rendererfg.SetActiveCamera(camera)
    renderWindow.Render()

    camera.OrthogonalizeViewUp();
    rendererfg.ResetCameraClippingRange();
    renderer.ResetCameraClippingRange();
    
    renderWindow.Render()
    windowToImageFilter = vtk.vtkWindowToImageFilter()
    #windowToImageFilter.SetInputBufferTypeToRGBA()
    windowToImageFilter.SetInput(renderWindow)
    windowToImageFilter.SetScale(1)
    windowToImageFilter.Update()

    writer = vtk.vtkPNGWriter()
    writer.SetFileName(output)
    writer.SetInputConnection(windowToImageFilter.GetOutputPort())
    writer.Write()



def screencapMulticlass(ifile, output, casenum="137"):
    '''from https://gist.github.com/pangyuteng/facd430d0d9761fc67fff4ff2e5fffc3 
    Save screencap of imagenumpy in path to be used afterwards
    '''

    reader = vtk.vtkXMLImageDataReader()
    reader.SetFileName(ifile)
    reader.Update()
    # print(reader.GetOutput().GetDimensions())
    #print(len(np.transpose(image.nonzero())))
    # colormapping = {"2": (0.74,0.08,0.08),
    #                 "1": (0.9,0.69,0.0),
    #                 "3": (0.0,0.45,0.5)}
    
        # colormapping = {"2": (0.15,0.32,0.45),
        #         "1": (0.9,0.8,0.7),
        #         "3": (1.0,0.4,0.15)}



    colormapping = {"2": (0.10,0.3,0.65),
                "1": (0.9,0.8,0.7),
                "3": (1.0,0.4,0.15),
                "4": (1.0,1.0,0.20)}
    opacitymapping = {"2": 0.8,
            "1": 0.8,
            "3": 0.8,
            "4":0.99}


    renderer = vtk.vtkRenderer()
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.SetOffScreenRendering(1)
    renderWindow.AddRenderer(renderer)
    renderWindow.SetAlphaBitPlanes(1)
    renderWindow.SetSize(2000,1500)
    imgrange = reader.GetOutput().GetScalarRange()
    print(imgrange)
    actors = []
    center = (0,0,0)
    for u in range(int(imgrange[0])+1, int(imgrange[1])+1):#[1,2,3,4]:
        threshold = vtk.vtkImageThreshold()
        threshold.SetInputConnection(reader.GetOutputPort()) #Connection(reader.GetOutputPort())
        threshold.ThresholdBetween(u, u)  #th
        threshold.ReplaceInOn()
        threshold.SetInValue(1)  # set all values below th to 0
        threshold.ReplaceOutOn()
        threshold.SetOutValue(0)  # set all values above th to 1
        threshold.Update()

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

        actor.GetProperty().SetColor(colormapping[str(u)])
        actor.GetProperty().SetOpacity(opacitymapping[str(u)])
        #actor.RotateX(90)

        renderer.AddActor(actor)
        center = actor.GetCenter()
    
    renderer.SetBackground(1.0, 1.0, 1.0)

    camera = renderer.MakeCamera()
    #camera.SetPosition(0,0,800)
    #camera.SetPosition(400,150,0) # spacing 0.227

    if casenum == "137":
        camera.SetPosition(1000,800,530) # spacing 1
        camera.SetFocalPoint(231.72304458265867, 281.3735192306137, 282.29167910976366) 
        camera.SetViewUp(-0.1, -0.28, 0.95)
    elif casenum == "138":
        camera.SetPosition(1220.5390876788997, 673.1567789943556, 335.78085552031774) # spacing 1
        camera.SetFocalPoint(249.32250085545363, 267.5176803725833, 261.3092929495787) 
        camera.SetViewUp(-0.0690215000169431, -0.03479958081559921, 0.9970080349277382)
    elif casenum == "139":
        camera.SetPosition(1234.205896352638, -117.83510206877108, 387.1811305239322) # spacing 1
        camera.SetFocalPoint(255.50000000000026, 255.50000000000009, 254.5000000000001) 
        camera.SetViewUp(-0.16507953959775667, -0.038287930104191355, 0.9855367978997687)
    elif casenum == "142":
        camera.SetPosition(1409.621676054513, 388.3356818559007, 357.9944059869261) # spacing 1
        camera.SetFocalPoint(256.1453510400097, 255.56808239880462, 247.6706778295724) 
        camera.SetViewUp(-0.09407334171245772, -0.009924426196404655, 0.9955158020562624)
    elif casenum == "143":
        camera.SetPosition(-650.7999857855324, 316.78600067829495, 155.43119583499293) # spacing 1
        camera.SetFocalPoint(254.80673971036933, 258.0625813138862, 220.46727595671342) 
        camera.SetViewUp(-0.07377299088327371, -0.03385072040657621, 0.9967003935707522)
    elif casenum == "145":
        camera.SetPosition(1272.3814126235534, -28.889241733865617, 341.44059747993856) # spacing 1
        camera.SetFocalPoint(255.50000000000023, 255.50000000000006, 254.5000000000002) 
        camera.SetViewUp(-0.12426108701663868, -0.0987889987325104, 0.9873195612277068)
    elif casenum == "146":
        camera.SetPosition(-480.64639232126416, 649.9945753572392, 217.0695093097312) # spacing 1
        camera.SetFocalPoint(261.25593485674057, 266.73312300993416, 132.0895318529527) 
        camera.SetViewUp(0.08499984246320676, -0.05603703493908549, 0.9948039392244412)
    elif casenum == "148":
        camera.SetPosition(-363.7641503684439, 1841.740680994073, 142.223206886471) # spacing 1
        camera.SetFocalPoint(255.5000000000001, 255.4999999999997, 253.99999999999974) 
        camera.SetViewUp(-0.1717149865567503, 0.0023825377747678234, 0.9851437899644736)
    elif casenum == "162":
        camera.SetPosition(912.1360403404894, 1219.4166282059773, 254.6929328224841) # spacing 1
        camera.SetFocalPoint(255.50000000000006, 255.49999999999997, 254.49999999999994) 
        camera.SetViewUp(-0.03700687572629638, 0.025009743955626475, 0.9990020039300485)     
    elif casenum == "190":
        camera.SetPosition(81.44595323522267, 1870.5127465854487, -272.1651215843878) # spacing 1
        camera.SetFocalPoint(255.49999999999986, 255.49999999999986, 254.50000000000023) 
        camera.SetViewUp(-0.07570822879169399, 0.3017627113640458, 0.9503723113198554)     
    elif casenum == "192":
        camera.SetPosition(136.61874409891658, 1358.5754365039545, -105.21936451361597) # spacing 1
        camera.SetFocalPoint(255.5, 255.5, 254.5) 
        camera.SetViewUp(-0.07570822879169399, 0.3017627113640458, 0.9503723113198554)    
    elif casenum == "163":
        # camera.SetPosition(991.2817211421429, -605.5551085868585, -23.91614830800691) # spacing 1
        # camera.SetFocalPoint(255.5000000000001, 255.49999999999957, 254.5000000000002) 
        # camera.SetViewUp(0.044394713655192354, -0.2728106929000861, 0.9610428893857176)
        camera.SetPosition(1082.843131000895, -152.90436517958975, -459.53254983613095) # spacing 1
        camera.SetFocalPoint(272.120231576054, 248.165519921031, 276.7944973668337) 
        camera.SetViewUp(0.2856371298967383, -0.6737763549117353, 0.6814960407707511)      
    elif casenum == "177":
        # camera.SetPosition(690.9347955520454, 1112.0667517889565, 626.2356060635063) # spacing 1
        # camera.SetFocalPoint(220.6152792316476, 312.97492897695184, 294.4309474200802) 
        # camera.SetViewUp(0.057999847002381484, -0.41176253845691324, 0.9094435824564487) 
        camera.SetPosition(330.06061945450085, -766.0283677763737, 153.57756716278692) # spacing 1
        camera.SetFocalPoint(223.52710984992112, 205.5976718635301, 273.7273513463957) 
        camera.SetViewUp(0.17620090464872576, -0.10175219057806838, 0.9790810655474527)      
    elif casenum == "194":
        # camera.SetPosition(690.9347955520454, 1112.0667517889565, 626.2356060635063) # spacing 1
        # camera.SetFocalPoint(220.6152792316476, 312.97492897695184, 294.4309474200802) 
        # camera.SetViewUp(0.057999847002381484, -0.41176253845691324, 0.9094435824564487) 
        # camera.SetPosition(1016.3155743032796, 874.2863804281911, 351.6890125632976) # spacing 1
        # camera.SetFocalPoint(235.5091650595668, 279.1853777201952, 273.98014515282676) 
        # camera.SetViewUp(0.0042712830451648726, -0.13498874609218414, 0.9908379254800495)  
        # camera.SetPosition(176.89317125378793, 119.44684296667383, 338.50862258646816) # spacing 1
        # camera.SetFocalPoint(269.0349198391431, 253.34569168733125, 268.11962547264505) 
        # camera.SetViewUp(-0.2819577714595236, 0.5912390800225928, 0.7556031798289478)  
        # camera.SetPosition(-928.0086345059518, 120.61989189153255, 287.2433403515401) # spacing 1
        # camera.SetFocalPoint(255.26266497455592, 256.7569724646287, 251.43606504964842) 
        # camera.SetViewUp(0.04653638264257264, -0.1445083474366278, 0.9884086718618362)      
        camera.SetPosition(-0.07771454318191218, 211.83170928157543, 272.90595920292253) # spacing 1
        camera.SetFocalPoint(257.436595768047, 241.4591044756921, 265.11325282397195) 
        camera.SetViewUp(0.04653638264257264, -0.1445083474366278, 0.9884086718618362)                                    

    renderer.SetActiveCamera(camera)
    renderWindow.Render()

    # print(position,focal_point,view_up,)

    # print(camera.GetViewTransformMatrix())
    # print(camera.GetViewTransformMatrix().GetElement(0,0))
    # print(camera.GetViewTransformMatrix().GetElement(0,1))
    # print(camera.GetViewTransformMatrix().GetElement(0,2))

    frame_list = []
    #for n,q in enumerate([60]*6):
    # n = 0
    # q = 60
    # transform = vtk.vtkTransform()
    # transform.Identity()

    # transform.Translate(*center)
    # transform.RotateWXYZ(q,view_up)
    # transform.RotateWXYZ(0,axis)
    # transform.Translate(*[-1*x for x in center])

    # new_position = [0,0,0]
    # new_focal_point = [0,0,0]
    # transform.TransformPoint(position,new_position)
    # transform.TransformPoint(focal_point,new_focal_point)

    # camera.SetPosition(new_position)
    # camera.SetFocalPoint(new_focal_point)

    # focal_point = camera.GetFocalPoint()
    # view_up = camera.GetViewUp()
    # position = camera.GetPosition() 

    camera.OrthogonalizeViewUp();
    renderer.ResetCameraClippingRange();
    
    renderWindow.Render()
    windowToImageFilter = vtk.vtkWindowToImageFilter()
    #windowToImageFilter.SetInputBufferTypeToRGBA()
    windowToImageFilter.SetInput(renderWindow)
    windowToImageFilter.SetScale(1)
    windowToImageFilter.Update()

    writer = vtk.vtkPNGWriter()
    writer.SetFileName(output)
    writer.SetInputConnection(windowToImageFilter.GetOutputPort())
    writer.Write()




def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-ipath", help="path to file or folder to capture", default="", required=True)
    parser.add_argument("-opath", help="path to folder to save screencap", default="", required=False)
    parser.add_argument("--all", help="capture all files in folder", action="store_true", required=False, default=False)
    parser.add_argument("-case", help="case", default="", required=False)


    args = parser.parse_args()
    ipath = args.ipath
    opath = args.opath
    allfiles = args.all
    case = args.case

    # if allfiles and isdir(ipath):
    #     thefiles = getsubfiles(ipath, suffix=".nii.gz")
    # elif isfile(ipath):
    #     thefiles = [ipath]
    # else:
    #     print("Unable to read files")
    #     return
    # team = ipath.split(os.path.sep)[-2]
    # print(team)
    # for f in thefiles:
    #     filename = team + "_" + (f.split(os.path.sep)[-1]).split(".")[0]
    #     print(filename)
    #     print(f)
    #     print(opath)
    #     screencap2(f, opath, filename)

    thefiles = getsubfiles(ipath, contains=case, suffix="web.vti")
    for f in thefiles:
        casenum = ((f.split(os.path.sep)[-1]).split(".")[0]).split("_")[0]
        print(casenum)
        #output = "/media/camila/Datos4TB/challenge/submissions/test-1/screencap/"+casenum+"_"+ipath.split(os.path.sep)[-1]+".png"
        output = "/media/camila/Datos4TB/challenge/submissions/test-2/screencap/"+casenum+"_"+ipath.split(os.path.sep)[-2]+"_f.png"
        print(output)
        screencapMulticlass(f, output, casenum)

if __name__ == "__main__":
    main()
