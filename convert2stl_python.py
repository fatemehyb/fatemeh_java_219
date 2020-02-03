# https://www.researchgate.net/post/How_do_I_convert_a_mesh_in_a_voxelized_volume_3d_image2
import SimpleITK as sitk
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import ConvexHull
import numpy as np
# import scipy.io as sio
# from stl import mesh
# # import shapely
# import random
# import mayavi
from mayavi import mlab
def ismember(A, B):
    return [ np.sum(a == B) for a in A ]


myvolume=np.random.randint(1,101, size=(30,3))
points=myvolume
w1 = myvolume[2]
h1 = myvolume[1]
d1 = myvolume[0]
# The_Normal,p1,p2,p3,index_changed=main_normal(myvolume,spacing,verts2,faces2)
hull=ConvexHull(myvolume)
# hull = ConvexHull(pts)

# fig = plt.figure()
# ax = fig.add_subplot(111, projection="3d")
#
# # Plot defining corner points
# ax.plot(((np.asarray(myvolume).transpose()))[0], ((np.asarray(myvolume).transpose()))[1], ((np.asarray(myvolume).transpose()))[2], "ko")
#
# # 12 = 2 * 6 faces are the simplices (2 simplices per square face)
# for s in hull.simplices:
#     s = np.append(s, s[0])  # Here we cycle back to the first coordinate
#     ax.plot(myvolume[s[0:3][0]], myvolume[s[0:3][1]], myvolume[s[0:3][2]], "r-")
#
# # Make axis label
# # for i in ["x", "y", "z"]:
# #     ["convex hull of points"][1]eval("ax.set_{:s}label('{:s}')".format(i, i))
# # plt.fill((hull.points.transpose())[0],(hull.points.transpose())[2],(hull.points.transpose())[1])
# plt.show()
# del(myvolume)
# hfacets=hull.simplices.shape[1]
# d=points.shape[0]
#
# for i in range(1,hfacets):
#     for j in range(i+1,hfacets):
#         if np.count_nonzero(ismember(hull.simplices[i,:],hull.simplices[j,:]))==d:
#             print('Duplicate hull facet.')

from matplotlib.tri import Triangulation
# tri = Triangulation(np.ravel(w), np.ravel(theta))

# ax = plt.axes(projection='3d')
from stl import mesh
import stl
mlab.triangular_mesh(points[:,0],points[:,1],points[:,2], hull.simplices)
# mesh.save('mesh.stl')
mlab.show()
cube = mesh.Mesh(np.zeros(hull.simplices.shape[0], dtype=mesh.Mesh.dtype))
for i, f in enumerate(hull.simplices):
    for j in range(3):
        cube.vectors[i][j] = hull.points[f[j],:]

# Write the mesh to file "cube.stl"
cube.save('cube2.stl')
############################running matlab scripts in python #################################
from oct2py import octave as oct
# octave = oct.oct2py('D:\Octave\Octave-5.1.0.0\mingw64\\bin\octave-cli.exe')
import os

oct.eval("cd D:\matlab_useful_codes\Mesh_voxelisation\Mesh_voxelisation")
cwd = os.getcwd()

oct.addpath(cwd)
oct.addpath('D:\matlab_useful_codes\Mesh_voxelisation\Mesh_voxelisation')
oct.eval("VOXELISE_example")
oct.eval("cd D:\matlab_useful_codes\Mesh_voxelisation")
oct.eval("save -v7 myworkspace.mat")

from scipy.io import loadmat
D = loadmat("D:\matlab_useful_codes\Mesh_voxelisation\myworkspace.mat")
print(D.keys())

########reading the .mat matrix convert it to numpy array and save it as .nrrd image using SimpleITK
# z=sio.loadmat('test_voxel.mat')
z2=D['OUTPUTgrid']
z2_Image=sitk.GetImageFromArray(z2)
sitk.WriteImage(z2_Image,'output.nrrd')
