import numpy as np
import sys

###################################################################################################
# This python script reads in an ascii file generated by the ASPECT gravity postprocessor and
# generates a vtu file of the measuring points with the gravity vector (gx,gy,gz) and potential U.
# 
# How to use it:
#
# > python3 gravity_ascii_to_vtu.py file_name
#
# If no file_name is provided this scripts uses gravity-00000 by default.
# The generated vtu file is named: file_name_sphere.vtu
#
# For questions/remarks: c.thieulot@uu.nl
###################################################################################################

if int(len(sys.argv) == 2):
   filename = sys.argv[1]
else:
   filename = 'gravity-00000'
   
print('processing file:',filename)

###################################################################################################
# automatic detection of number of points
###################################################################################################

f = open(filename, 'r')
lines = f.readlines()
f.close
counter=0

counter=0
for i in lines:
    counter+=1

nheader=30 # number of header lines in gravity file

nnp=counter-nheader

print('file counts',counter,'lines')
print('header counts',nheader,'lines')
print('number of data points is then',nnp)

###################################################################################################
# read in values from file
###################################################################################################

r  = np.zeros(nnp,dtype=np.float64)  
t  = np.zeros(nnp,dtype=np.float64)  
p  = np.zeros(nnp,dtype=np.float64)  
x  = np.zeros(nnp,dtype=np.float64)  
y  = np.zeros(nnp,dtype=np.float64)  
z  = np.zeros(nnp,dtype=np.float64)  
gx = np.zeros(nnp,dtype=np.float64)  
gy = np.zeros(nnp,dtype=np.float64)  
gz = np.zeros(nnp,dtype=np.float64)  
U  = np.zeros(nnp,dtype=np.float64)  

for i in range(0,nnp):
       vals=lines[i+nheader].strip().split()
       r[i]=vals[0]  
       t[i]=vals[2]  
       p[i]=vals[1]  
       x[i]=vals[3]  
       y[i]=vals[4]  
       z[i]=vals[5]  
       gx[i]=vals[6]  
       gy[i]=vals[7]  
       gz[i]=vals[8]  
       U[i]=vals[11]  

###################################################################################################

vtufile=open(filename+"_sphere.vtu","w")
vtufile.write("<VTKFile type='UnstructuredGrid' version='0.1' byte_order='BigEndian'> \n")
vtufile.write("<UnstructuredGrid> \n")
vtufile.write("<Piece NumberOfPoints=' %5d ' NumberOfCells=' %5d '> \n" %(nnp,nnp))
#####
vtufile.write("<Points> \n")
vtufile.write("<DataArray type='Float32' NumberOfComponents='3' Format='ascii'> \n")
for i in range(0,nnp):
       vtufile.write("%10f %10f %10f \n" %(x[i],y[i],z[i]))
vtufile.write("</DataArray>\n")
vtufile.write("</Points> \n")
#####
vtufile.write("<PointData>\n")
vtufile.write("<DataArray type='Float32' NumberOfComponents='3' Name='gravity vector g' Format='ascii'> \n")
for i in range(0,nnp):
       vtufile.write("%10e %10e %10e \n" %(gx[i],gy[i],gz[i]))
vtufile.write("</DataArray>\n")
vtufile.write("<DataArray type='Float32' Name='r' Format='ascii'> \n")
for i in range(0,nnp):
       vtufile.write("%10e \n" %(r[i]))
vtufile.write("</DataArray>\n")
vtufile.write("<DataArray type='Float32' Name='theta' Format='ascii'> \n")
for i in range(0,nnp):
       vtufile.write("%10e \n" %(t[i]))
vtufile.write("</DataArray>\n")
vtufile.write("<DataArray type='Float32' Name='phi' Format='ascii'> \n")
for i in range(0,nnp):
       vtufile.write("%10e \n" %(p[i]))
vtufile.write("</DataArray>\n")
vtufile.write("<DataArray type='Float32' Name='U' Format='ascii'> \n")
for i in range(0,nnp):
       vtufile.write("%10e \n" %(U[i]))
vtufile.write("</DataArray>\n")
#--  
vtufile.write("</PointData>\n")
#####
vtufile.write("<Cells>\n")
vtufile.write("<DataArray type='Int32' Name='connectivity' Format='ascii'> \n")
for i in range (0,nnp):
       vtufile.write("%d \n" % i) 
vtufile.write("</DataArray>\n")
vtufile.write("<DataArray type='Int32' Name='offsets' Format='ascii'> \n")
for i in range (0,nnp):
       vtufile.write("%d \n" %(i+1))
vtufile.write("</DataArray>\n")
vtufile.write("<DataArray type='Int32' Name='types' Format='ascii'>\n")
for iel in range (0,nnp):
       vtufile.write("%d \n" % 1) 
vtufile.write("</DataArray>\n")
vtufile.write("</Cells>\n")
#####
vtufile.write("</Piece>\n")
vtufile.write("</UnstructuredGrid>\n")
vtufile.write("</VTKFile>\n")
vtufile.close()

print("generated", filename+'_sphere.vtu')

###################################################################################################