import read_file_function as read
import sys,os
import SimpleITK as sitk
import nibabel as nib
string_write="D:/fatemeh/UWO_CASES/UWO_CASE_1649R/1649R_154um_NRRD_Segmentations/1649R_154um-label-SS.nrrd"

# os.path.split(string_write)[0]
x=string_write.split("/")
img=sitk.ReadImage(string_write)
size = img.GetSize()
# img=read.read(string_write)
array1=sitk.GetArrayFromImage(img)
img2=sitk.GetImageFromArray(array1)
sitk.WriteImage(img2,"D:/Fatemeh/UWO_CASES/Labels"+"//"+x[-1]+".nii.gz")
