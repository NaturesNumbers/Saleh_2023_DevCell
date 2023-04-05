#!/usr/bin/env python
# -*- coding: latin-1 -*-
# Based on a code by Emmanuel Faure

import sys,os
from skimage.io import imsave
import vtk
from vtk.util.numpy_support import vtk_to_numpy
from skimage import io
import numpy as np
from os import listdir
from os.path import isfile, join

#Convert in OBJ
def getObjSurface(ndata,border=4,voxel=[1,1,1]): #NDATA is boolean 
	pixels= np.array(np.where(ndata))
	xmin=pixels[0,:].min()
	xmax=pixels[0,:].max()
	ymin=pixels[1,:].min()
	ymax=pixels[1,:].max()
	zmin=pixels[2,:].min()
	zmax=pixels[2,:].max()

	pixels[0,:]-=xmin
	pixels[1,:]-=ymin
	pixels[2,:]-=zmin

	data=np.zeros([2*border+xmax-xmin,2*border+ymax-ymin,2*border+zmax-zmin],dtype=np.uint8)
	data[pixels[0,:]+border,pixels[1,:]+border,pixels[2,:]+border]=255

	vtk_mesh = vtk.vtkPolyData()
	vtk_points = vtk.vtkPoints()
	vtk_triangles = vtk.vtkCellArray()
	vtk_cells = vtk.vtkLongArray()

	data_string = data.tobytes('F')

	reader = vtk.vtkImageImport()
	reader.CopyImportVoidPointer(data_string, len(data_string))
	reader.SetDataScalarTypeToUnsignedChar()
	reader.SetNumberOfScalarComponents(1)
	reader.SetDataSpacing(voxel)
	reader.SetDataExtent(0, data.shape[0] - 1, 0, data.shape[1] - 1, 0, data.shape[2] - 1)
	reader.SetWholeExtent(0, data.shape[0] - 1, 0, data.shape[1] - 1, 0, data.shape[2] - 1)
	reader.Update()

	image=reader.GetOutput()
	image.SetSpacing(voxel)

	#MARCHING CUBES
	contour = vtk.vtkDiscreteMarchingCubes()
	contour.SetInputData(image)
	contour.ComputeNormalsOn()
	contour.ComputeGradientsOn()
	contour.SetValue(0,255)
	contour.Update()

	polydata= contour.GetOutput()

	obj=""
	if not polydata.GetPoints() is None :
	    for p in range(polydata.GetPoints().GetNumberOfPoints()):
	        v=polydata.GetPoints().GetPoint(p) 
	        obj+='v ' + str((v[0]+xmin-border)*cratio) +' '+str(v[1]+ymin-border) +' '+str(v[2]+zmin-border)+'\n'
	    for t in range(polydata.GetNumberOfCells()):
	         obj+='f ' + str(1+polydata.GetCell(t).GetPointIds().GetId(0)) +' '+str(1+polydata.GetCell(t).GetPointIds().GetId(1)) +' '+str(1+polydata.GetCell(t).GetPointIds().GetId(2))+'\n'
	return obj
 

def saveOBJ(obj,filename):
    f=open(filename,'w')
    f.write(obj)
    f.close()

			
########################## MAIN 

mypath = input("What is the path to your masks? ")
cratio = float(input("What is the slice to pixel ratio? "))
onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f)) and f[-3::]=='tif' ]

for f in onlyfiles:
	cell_matrix = io.imread(mypath+'/'+f)
	ndata=np.zeros(cell_matrix.shape,np.bool8)
	ndata=cell_matrix
	obj=getObjSurface(ndata)
	saveOBJ(obj,mypath+'/'+f[0:-3]+'obj')