#version date: 8/11/2018

import fileinput
import math
import numpy as np
import sys
import os

sp="    "
xcells=150
ycells=xcells
zcells=115
xsize=15
ysize=xsize
dx=xsize*2/xcells
dy=dx

Btype=0 #0: classical+chiral, 1: classical only, 2: chiral only
Blimit=1.e-6

if(not((Btype == 0) or (Btype == 1) or (Btype == 2))):
  print("B type option unclear... I quit.\n")
  sys.exit(1)

ix=1-1
iy=2-1
iz=3-1
iH1x=4-1
iH2x=6-1
iH1Cx=8-1
iH2Cx=10-1
iH1y=12-1
iH2y=14-1
iH1Cy=16-1
iH2Cy=18-1
iH1Cz=20-1
iH2Cz=22-1

xpoints_all=np.linspace(-xsize+dx/2,xsize-dx/2,num=xcells,endpoint=True)
ypoints_all=np.linspace(-ysize+dy/2,ysize-dy/2,num=ycells,endpoint=True)
xpoints=np.round(xpoints_all[::5],decimals=2)
ypoints=np.round(ypoints_all[::5],decimals=2)
zcoord=0

if(len(sys.argv)<3):
   print ('Syntax: ./prepare_arrows.py <inputfile> <outputfile>')
   sys.exit(2)

inputfile=sys.argv[1]
outputfile=sys.argv[2]

datein = open(inputfile,"r")
dateo = open(outputfile,"w")

dateo.write("#1 x #2 y #3 H1x #4 H1y #5 H1 #6 H2x #7 H2y #8 H2 #9 H1Cx #10 H1Cy #11 H1C #12 H2Cx #13 H2Cy #14 H2C #15 Hxt #16 Hyt #17 Ht\n")

for line in datein:
  stuff=line.split()
  if(len(stuff)>0):
   if(stuff[0][0] != "#"):
    z=float(stuff[iz])
    if(z == zcoord):
     x=np.round(float(stuff[ix]),decimals=2)
     if x in xpoints[:]:
      y=np.round(float(stuff[iy]),decimals=2)
      if y in ypoints[:]:
         H1x_d=float(stuff[iH1x])
         H2x_d=float(stuff[iH2x])
         H1Cx_d=float(stuff[iH1Cx])
         H2Cx_d=float(stuff[iH2Cx])
         H1y_d=float(stuff[iH1y])
         H2y_d=float(stuff[iH2y])
         H1Cy_d=float(stuff[iH1Cy])
         H2Cy_d=float(stuff[iH2Cy])
         H1Cz_d=float(stuff[iH1Cz])
         H2Cz_d=float(stuff[iH2Cz])

         H1m=math.sqrt(H1x_d**2+H1y_d**2)
         if(H1m>Blimit):
           H1x=H1x_d/H1m
           H1y=H1y_d/H1m
         else:
           H1x=0
           H1y=0
           H1m=0
    
         H2m=math.sqrt(H2x_d**2+H2y_d**2)
         if(H2m>Blimit):
           H2x=H2x_d/H2m
           H2y=H2y_d/H2m
         else:
           H2x=0
           H2y=0
           H2m=0

         H1Cm=math.sqrt(H1Cx_d**2+H1Cy_d**2)
         if(H1Cm>Blimit):
           H1Cx=H1Cx_d/H1Cm
           H1Cy=H1Cy_d/H1Cm
         else:
           H1Cx=0
           H2Cx=0
           H1Cm=0

         H2Cm=math.sqrt(H2Cx_d**2+H2Cy_d**2)
         if(H2Cm>Blimit):
           H2Cx=H2Cx_d/H2Cm
           H2Cy=H2Cy_d/H2Cm
         else:
           HC2x=0
           HC2y=0
           HC2m=0

         Hxt=H1x_d+H2x_d+H1Cx_d+H2Cx_d
         Hyt=H1y_d+H2y_d+H1Cy_d+H2Cy_d
         Hm=math.sqrt(Hxt**2+Hyt**2)
         if(Hm>Blimit):
           Hxt=Hxt/Hm
           Hyt=Hyt/Hm
         else:
           Hxt=0
           Hyt=0

         dateo.write('{:5.2f}'.format(x)+sp+'{:5.2f}'.format(y)+sp+'{:14.8e}'.format(H1x)+sp+'{:14.8e}'.format(H1y)+sp+'{:14.8e}'.format(H1m)+sp+'{:14.8e}'.format(H2x)+sp+'{:14.8e}'.format(H2y)+sp+'{:14.8e}'.format(H2m)+sp+'{:14.8e}'.format(H1Cx)+sp+'{:14.8e}'.format(H1Cy)+sp+'{:14.8e}'.format(H1Cm)+sp+'{:14.8e}'.format(H2Cx)+sp+'{:14.8e}'.format(H2Cy)+sp+'{:14.8e}'.format(H2Cm)+sp+'{:14.8e}'.format(Hxt)+sp+'{:14.8e}'.format(Hyt)+sp+'{:14.8e}'.format(Hm)+"\n")
 

dateo.close()
print("All done.\n")
