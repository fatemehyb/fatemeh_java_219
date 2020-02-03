import read_file_function as read
import sys,os
import SimpleITK as sitk
import nibabel as nib
string_write="D:/fatemeh/UWO_CASES/UWO_CASE_1649R/1649R_154um_DICOM"

# os.path.split(string_write)[0]
x=string_write.split("/")
reader = sitk.ImageSeriesReader()

dicom_names = reader.GetGDCMSeriesFileNames( string_write )
reader.SetFileNames(dicom_names)

img = reader.Execute()

size = img.GetSize()
# img=read.read(string_write)
array1=sitk.GetArrayFromImage(img)
img2=sitk.GetImageFromArray(array1)
sitk.WriteImage(img2,"D:/Fatemeh/UWO_CASES"+"//"+x[-1]+".nii.gz")
