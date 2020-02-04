
#This program performs surface and volume rendering for a dataset.
import vtk
import itk
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
import argparse
import numpy as np

#Gaussian filter
def gaussianFilter(image):

    filter = vtk.vtkImageGaussianSmooth()

    filter.SetInputConnection(image.GetOutputPort())
    filter.SetStandardDeviation(1)
    filter.SetRadiusFactors(1,1,1)
    filter.SetDimensionality(3)
    filter.Update()

    return filter.GetOutput()

# Thresholding
def threshold(image, lowerThreshold, upperThreshold):

    thresh = vtk.vtkImageThreshold()
    thresh.SetInputData(image)

    thresh.ThresholdBetween(lowerThreshold, upperThreshold)
    thresh.ReplaceInOn()
    thresh.SetInValue(1)
    thresh.ReplaceOutOn()
    thresh.SetOutValue(0)
    thresh.Update()

    return thresh.GetOutput()

# Creates mesh using marching cubes for an image
def createMesh(mask, threshold1,threshold2):

    mesh = vtk.vtkMarchingCubes()
    mesh.SetInputData(mask)
    mesh.SetValue(threshold1, threshold2)
    mesh.Update()

    return mesh.GetOutput()

# Creates mesh for a thresholded image
def createThresholdMesh(mask,threshold1,threshold2):

    mesh = vtk.vtkDiscreteMarchingCubes()
    mesh.SetInputData(mask)
    mesh.SetValue(threshold1,threshold2)
    # mesh.GenerateValues(threshold2-threshold1+1,threshold1,threshold2)
    mesh.Update()

    return mesh.GetOutput()

def createThresholdMesh2(mask,threshold1,threshold2):

    mesh = vtk.vtkDiscreteMarchingCubes()
    mesh.SetInputData(mask)
    # mesh.SetValue(threshold1,threshold2)
    mesh.GenerateValues(threshold2-threshold1+1,threshold1,threshold2)
    mesh.Update()

    return mesh.GetOutput()



# Extracts two isosurfaces representing skin and bone
def dualSurfaceRendering(image1):

    ren = vtk.vtkRenderer()
    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(ren)
    renWin.PolygonSmoothingOn()
    renWin.SetSize(400,400)

    # create a renderwindowinteractor
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)

    # skinMesh = createMesh(image1,(int)(image1.GetScalarRange()[0]),(int)(image1.GetScalarRange()[1]-1))
    skinMesh = createMesh(image1,0,1)
    # boneMesh = createThresholdMesh(image2,0, 255)
    ##########################################################
    # # polydata=vtk.vtkPolyData()
    # # polydata.SetPoints(boneMesh.GetPoints())
    # splatter=vtk.vtkGaussianSplatter()
    # splatter.SetInputData(boneMesh)
    # splatter.SetRadius(0.02)
    # surface=vtk.vtkContourFilter()
    # surface.SetInputConnection(splatter.GetOutputPort())
    # surface.SetValue(0,0.01)
    # # # Convert the image to a polydata
    # # imageDataGeometryFilter = vtk.vtkImageDataGeometryFilter()
    # # imageDataGeometryFilter.SetInputData(image2)
    # # imageDataGeometryFilter.Update()
    # mapper = vtk.vtkPolyDataMapper()
    # mapper.SetInputConnection(surface.GetOutputPort())
    # actor = vtk.vtkActor()
    # actor.SetMapper(mapper)
    # # actor.GetProperty().SetPointSize(0.01)
    #############################################################

    # boneMesh = createThresholdMesh(image1,(int)(image1.GetScalarRange()[1]-1), (int)(image1.GetScalarRange()[1]))
    skinMapper = vtk.vtkPolyDataMapper()
    skinMapper.SetInputData(skinMesh)
    skinMapper.ScalarVisibilityOff()

    skinActor = vtk.vtkActor()
    skinActor.SetMapper(skinMapper)
    skinActor.GetProperty().SetColor(1.0,1.0,1.0)
    skinActor.GetProperty().SetOpacity(1.0)
    skinActor.GetProperty().SetAmbient(0.1)
    skinActor.GetProperty().SetDiffuse(0.9)
    skinActor.GetProperty().SetSpecular(0.1)

    # boneMapper = vtk.vtkPolyDataMapper()
    # boneMapper.SetInputData(boneMesh)
    # boneMapper.ScalarVisibilityOff()
    #
    # boneActor = vtk.vtkActor()
    # boneActor.SetMapper(boneMapper)
    # boneActor.GetProperty().SetColor(1.0,0.0,0.0)
    # boneActor.GetProperty().SetAmbient(0.2)
    # boneActor.GetProperty().SetDiffuse(0.7)
    # boneActor.GetProperty().SetSpecular(0.2)
    # boneActor.GetProperty().SetOpacity(1.0)
    # boneActor.GetProperty().SetPointSize(100)




    # sphere_view = vtk.vtkSphereSource()
    # # sphere_view.SetOrigin([-48.96,-48.96,-48.42])
    # # sphere_view.SetSpacing(0.18,0.18,0.18)
    # sphere_view.SetCenter(int(-72*0.18),int(441*0.18),int(439*0.18))
    # # sphere_view.SetCenter(int(250*0.18),int(150*0.18),int(600*0.18))
    # sphere_view.SetRadius(1.0)
    # # sphere_view.SetColor(0.0,1.0,1.0)
    # # mapper
    # mapper = vtk.vtkPolyDataMapper()
    # if vtk.VTK_MAJOR_VERSION <= 5:
    #    mapper.SetInput(sphere_view.GetOutput())
    # else:
    #    mapper.SetInputConnection(sphere_view.GetOutputPort())
    #
    # # actor
    # sphere_actor = vtk.vtkActor()
    # sphere_actor.SetMapper(mapper)
    # sphere_actor.GetProperty().SetColor(1.0,0.0,1.0)


    # ren.AddActor(boneActor)
    ren.AddActor(skinActor)
    # ren.AddActor(sphere_actor)
    # ren.AddActor(sigmoidActor)
    # ren.AddActor(actor)

    iren.Initialize()
    renWin.Render()
    iren.Start()
    stlWriter = vtk.vtkSTLWriter()
    stlWriter.SetFileName(str(("test2_30_october.stl")))
    stlWriter.SetInputData(skinMesh)
    stlWriter.Write()





# Reads, smooths and renders the dataset
def rendering(
        dataPath,
        dualSurface):

    import SimpleITK as sitk


    #
    # writing_path='data.nii'
    # reader = itk.ImageFileReader.New(FileName=dataPath)
    # reader.SetFileName(dataPath)
    # reader.Update()
    # image_input=reader.GetOutput()
    #
    # # imagee_out=sitk.GetImageFromArray(imagee)
    #
    # writer =itk.ImageFileWriter.New()
    # writer.SetFileName(writing_path)
    # writer.SetInput(image_input)
    # writer.Update()
    img=sitk.ReadImage(dataPath)
    img.SetSpacing([0.18,0.18,0.18])
    sitk.WriteImage(img,dataPath)
    # Read the data
    if dataPath.endswith('.nii'):
        reader = vtk.vtkNIFTIImageReader()
        # reader.SetSpacing(1,1,1)
        reader.SetFileName(dataPath)
        reader.Update()
    elif dataPath.endswith('.nhdr') or dataPath.endswith('.nrrd'):
        reader = vtk.vtkNrrdReader()

        reader.SetFileName(dataPath)
        # reader.SetSpacing(0.18,0.18,0.18)
        reader.Update()
    else:
        reader = vtk.vtkDICOMImageReader()
        reader.SetDirectoryName(dataPath)
        reader.Update()





    # Smooth the image
    filteredImage = gaussianFilter(reader)

    # Surface rendering


    if dualSurface:
        dualSurfaceRendering(filteredImage)




# Argument Parsing
# parser = argparse.ArgumentParser(
#         description = """Surface and Volume rendering""")

def visualize(path_string):
    # Data file path
    # parser.add_argument("dataPath")
    dataPath=path_string
    # Label=label
    # sigmoid=sigmoid_sinus

    lowerThreshold=0

    upperThreshold=255

    surface=0

    dualSurface=1

    volume=0
    dual_volume=0



    rendering(dataPath,
    dualSurface)

parser = argparse.ArgumentParser(
    description = """visualizing the results""")
parser.add_argument("-args0", type = str, default = (("\\\\samba.cs.ucalgary.ca\\fatemeh.yazdanbakhsh\Documents\medical_imaging\Results\label_3d_9_3_oct_2_test.nrrd")),
    help = "address of label in .nii format")
parser.add_argument("-args1", type = str, default = ((('//samba.cs.ucalgary.ca/fatemeh.yazdanbakhsh/Documents/Data_Sets/Calgary/TBone-2015/TBoneCBCT-2015-10segmentations/L2963Left/Bone.nrrd'))),
    help = "address of volume in .nii format")
parser.add_argument("-args2", type = int, default = 0,
    help = "lower threshold")
parser.add_argument("-args3", type = int, default = 255,
    help = "upper threshold")
parser.add_argument("-args4", type = int, default = 0,
    help = "surface rendering")
parser.add_argument("-args5", type = int, default = 1,
    help = "dual surface rendering")
parser.add_argument("-args6", type = int, default = 0,
    help = "volume rendering")
parser.add_argument("-args7", type = int, default = 0,
    help = "dual volume rendering")
parser.add_argument("-args8", type = str, default = (("\\\\samba.cs.ucalgary.ca\\fatemeh.yazdanbakhsh\Documents\medical_imaging\Results\sig.nii")),
    help = "sigmoid sinus directory")
args = parser.parse_args()
# label=args.args0
lowerThreshold=args.args2
upperThreshold=args.args3
# surface=args.args4
dualSurface=args.args5
# volume=args.args6
# dual_voulme=args.args7
# sigmoid_sinus=args.args8
visualize(args.args1)
# visualize("U:\Documents\Data_Sets\Calgary\TBone-2015\TBoneCBCT-2015-10\L2963_L_modified_1_oct_2018",label)
