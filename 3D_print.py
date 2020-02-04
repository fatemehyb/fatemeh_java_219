
#This program performs surface and volume rendering for a dataset.
###########import vtk
import SimpleITK as sitk
import itk
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
import argparse
import numpy as np
import bpy
Sspacing=[0.154,0.154,0.154]
#Gaussian filter
def remove_island():


    context = bpy.context
    # edit mode
    bpy.ops.object.mode_set(mode='EDIT')
    # split into loose parts
    bpy.ops.mesh.separate(type='LOOSE')
    # object mode
    bpy.ops.object.mode_set(mode='OBJECT')

    parts = context.selected_objects
    # sort by number of verts (last has most)
    parts.sort(key=lambda o: len(o.data.vertices))
    # print
    for part in parts:
        print(part.name, len(part.data.vertices))
    # pop off the last
    parts.pop()
    # remove the rest
    for o in parts:
        bpy.data.objects.remove(o)


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

def dualSurfaceRendering_two(image1,image2,stl_address):

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
    boneMesh = createMesh(image2,0,1)
    ##########################################################
    # polydata=vtk.vtkPolyData()
    # polydata.SetPoints(boneMesh.GetPoints())
    splatter=vtk.vtkGaussianSplatter()
    splatter.SetInputData(boneMesh)
    splatter.SetRadius(0.02)
    surface=vtk.vtkContourFilter()
    surface.SetInputConnection(splatter.GetOutputPort())
    surface.SetValue(0,0.01)
    # # Convert the image to a polydata
    # imageDataGeometryFilter = vtk.vtkImageDataGeometryFilter()
    # imageDataGeometryFilter.SetInputData(image2)
    # imageDataGeometryFilter.Update()
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(surface.GetOutputPort())
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    # actor.GetProperty().SetPointSize(0.01)
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

    boneMapper = vtk.vtkPolyDataMapper()
    boneMapper.SetInputData(boneMesh)
    boneMapper.ScalarVisibilityOff()

    boneActor = vtk.vtkActor()
    boneActor.SetMapper(boneMapper)
    boneActor.GetProperty().SetColor(1.0,0.0,0.0)
    boneActor.GetProperty().SetAmbient(0.2)
    boneActor.GetProperty().SetDiffuse(0.7)
    boneActor.GetProperty().SetSpecular(0.2)
    boneActor.GetProperty().SetOpacity(1.0)
    boneActor.GetProperty().SetPointSize(100)




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


    ren.AddActor(boneActor)
    ren.AddActor(skinActor)
    # ren.AddActor(sphere_actor)
    # ren.AddActor(sigmoidActor)
    # ren.AddActor(actor)

    iren.Initialize()
    renWin.Render()
    iren.Start()
    stlWriter = vtk.vtkSTLWriter()
    stlWriter.SetFileName(str((stl_address)))
    stlWriter.SetInputData(skinMesh)
    stlWriter.Write()





# Extracts two isosurfaces representing skin and bone
def dualSurfaceRendering(image1,stl_address):

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
    stlWriter.SetFileName(str((stl_address)))
    stlWriter.SetInputData(skinMesh)
    stlWriter.Write()





# Reads, smooths and renders the dataset
def rendering(
        dataPath,
        dualSurface,img1_path,img2_path,img3_path,img4_path):

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
    img_facial=sitk.ReadImage(img1_path)
    img_facial.SetSpacing(Sspacing)

    img_sigmoid=sitk.ReadImage(img2_path)
    img_sigmoid.SetSpacing(Sspacing)

    img_tensor=sitk.ReadImage(img3_path)
    img_tensor.SetSpacing(Sspacing)

    img_inner=sitk.ReadImage(img4_path)
    img_inner.SetSpacing(Sspacing)
    segmented_path="segmented_path.nrrd"
    img_segmented=sitk.AddImageFilter()
    img_1=img_segmented.Execute(img_facial,img_sigmoid)
    img_segmented2=sitk.AddImageFilter()
    img_2=img_segmented2.Execute(img_tensor,img_inner)
    img_segmented3=sitk.AddImageFilter()
    img_3=img_segmented3.Execute(img_1,img_2)
    img_3.SetOrigin([0,0,0])
    sitk.WriteImage(img_3,segmented_path)


    img=sitk.ReadImage(dataPath)
    img.SetSpacing(Sspacing)
    img_subtract=sitk.SubtractImageFilter()
    img=img_subtract.Execute(img,img_3)
    img.SetSpacing(Sspacing)
    img.SetOrigin([0,0,0])
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

    if segmented_path.endswith('.nii'):
        reader2 = vtk.vtkNIFTIImageReader()

        # reader.SetSpacing(1,1,1)
        reader2.SetFileName(segmented_path)
        reader2.Update()
    elif segmented_path.endswith('.nhdr') or dataPath.endswith('.nrrd'):
        reader2 = vtk.vtkNrrdReader()

        reader2.SetFileName(segmented_path)
        # reader.SetSpacing(0.18,0.18,0.18)
        reader2.Update()
    else:
        reader2 = vtk.vtkDICOMImageReader()
        reader2.SetDirectoryName(segmented_path)
        reader2.Update()



    stl1="part1_L1537.stl"
    stl2="part2_L1537.stl"

    # Smooth the image
    filteredImage = gaussianFilter(reader)

    # Surface rendering


    if dualSurface:
        dualSurfaceRendering(filteredImage,stl1)


    # Smooth the image
    filteredImage2 = gaussianFilter(reader2)

    # Surface rendering


    if dualSurface:
        dualSurfaceRendering_two(filteredImage,filteredImage2,stl2)





# Argument Parsing
# parser = argparse.ArgumentParser(
#         description = """Surface and Volume rendering""")

def visualize(path_string,path1,path2,path3,path4):
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
    dualSurface,path1,path2,path3,path4)

parser = argparse.ArgumentParser(
    description = """visualizing the results""")
# parser.add_argument("-args0", type = str, default = (("//samba.cs.ucalgary.ca/fatemeh.yazdanbakhsh/Documents/medical_imaging/Results/label_3d_9_3_oct_2_test.nrrd")),
#     help = "address of label in .nii format")
# parser.add_argument("-args1", type = str, default = ((('//samba.cs.ucalgary.ca/fatemeh.yazdanbakhsh/Documents/Data_Sets/Calgary/TBone-2015/TBoneCBCT-2015-10segmentations/L2963Left/Bone.nrrd'))),
#     help = "address of volume in .nii format")
# parser.add_argument("-args9", type = str, default = ((('//samba.cs.ucalgary.ca/fatemeh.yazdanbakhsh/Documents/Data_Sets/Calgary/TBone-2015/TBoneCBCT-2015-10segmentations/L2963Left/FacialNerve.nrrd'))),
#     help = "address of volume in .nii format")
# parser.add_argument("-args10", type = str, default = ((('//samba.cs.ucalgary.ca/fatemeh.yazdanbakhsh/Documents/Data_Sets/Calgary/TBone-2015/TBoneCBCT-2015-10segmentations/L2963Left/SigmoidSinus.nrrd'))),
#     help = "address of volume in .nii format")
# parser.add_argument("-args11", type = str, default = ((('//samba.cs.ucalgary.ca/fatemeh.yazdanbakhsh/Documents/Data_Sets/Calgary/TBone-2015/TBoneCBCT-2015-10segmentations/L2963Left/TensorTympani.nrrd'))),
#     help = "address of volume in .nii format")
# parser.add_argument("-args12", type = str, default = ((('//samba.cs.ucalgary.ca/fatemeh.yazdanbakhsh/Documents/Data_Sets/Calgary/TBone-2015/TBoneCBCT-2015-10segmentations/L2963Left/InnerEar_L2963.nrrd'))),
#     help = "address of volume in .nii format")
# parser.add_argument("-args0", type = str, default = (("//samba.cs.ucalgary.ca/fatemeh.yazdanbakhsh/Documents/medical_imaging/Results/label_3d_9_3_oct_2_test.nrrd")),
#     help = "address of label in .nii format")
parser.add_argument("-args1", type = str, default = ((('D:/fatemeh/UWO_CASES/UWO_CASE_1537L/thereshould_L1537.nrrd'))),
    help = "address of volume in .nii format")
parser.add_argument("-args9", type = str, default = ((('D:/fatemeh/UWO_CASES/UWO_CASE_1537L/1537L_154um_NRRD_Segmentations/1537L_154um-label-FN.nrrd'))),
    help = "address of volume in .nii format")
parser.add_argument("-args10", type = str, default = ((('D:/fatemeh/UWO_CASES/UWO_CASE_1537L/1537L_154um_NRRD_Segmentations/1537L_154um-label-SS.nrrd'))),
    help = "address of volume in .nii format")
parser.add_argument("-args11", type = str, default = ((('D:/fatemeh/UWO_CASES/UWO_CASE_1537L/1537L_154um_NRRD_Segmentations/1537L_154um-label-ICA.nrrd'))),
    help = "address of volume in .nii format")
parser.add_argument("-args12", type = str, default = ((('D:/fatemeh/UWO_CASES/UWO_CASE_1537L/1537L_154um_NRRD_Segmentations/1537L_154um-label-IE.nrrd'))),
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
visualize(args.args1,args.args9,args.args10,args.args11,args.args12)
# visualize("U:\Documents\Data_Sets\Calgary\TBone-2015\TBoneCBCT-2015-10\L2963_L_modified_1_oct_2018",label)
