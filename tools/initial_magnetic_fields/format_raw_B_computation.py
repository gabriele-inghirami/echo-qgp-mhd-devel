# Author: Gabriele Inghirami (g.inghirami@gsi.de)
# License: PUBLIC DOMAIN

import fileinput
import math
import numpy as np
import sys
import os

sp="    "

Btype=0 #0: classical+chiral, 1: classical only, 2: chiral only
Blimit=1.e-6
convfactor=4.8113 #to convert from GeV^1/2/fm^3/2 to units of neutral pion mass squared

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

if(len(sys.argv)<3):
   print ('Syntax: ./extract.py <inputfile> <outputfile>')
   sys.exit(2)

inputfile=sys.argv[1]
outputfile=sys.argv[2]
outputfile4gp=sys.argv[2]+"_4gp"

datein = open(inputfile,"r")
dateoe = open(outputfile,"w")
dateog = open(outputfile4gp,"w")

dateog.write("#1 x #2 y #3 z #4 H1x #5 H2x #6 Hx #7 H1Cx #8 H2Cx #9 HCx #10 Hx+HCx #11 H1y #12 H2y #13 Hy #14 H1Cy #15 H2Cy #16 HCy #17 Hy+Hcy #18 H1Cz #19 H2Cz #20 HCz\n")

for line in datein:
  stuff=line.split()
  if(len(stuff)>0):
   if(stuff[0][0] != "#"):
    x=stuff[ix]
    y=stuff[iy]
    z=stuff[iz]
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
    
    H1x=H1x_d if abs(H1x_d) > Blimit else 0
    H2x=H2x_d if abs(H2x_d) > Blimit else 0
    Hx_d=H1x_d+H2x_d
    Hx=Hx_d if abs(Hx_d) > Blimit else 0
    H1Cx=H1Cx_d if abs(H1Cx_d) > Blimit else 0
    H2Cx=H2Cx_d if abs(H2Cx_d) > Blimit else 0
    HCx_d=H1Cx_d+H2Cx_d
    HCx=HCx_d if abs(HCx_d) > Blimit else 0
    if(Btype==0):
      Hxtot_d=H1x_d+H2x_d+H1Cx_d+H2Cx_d
    elif(Btype==1):
      Hxtot_d=Hx
    elif(Btype==2):
      Hxtot_d=Hxc
      
    Hxtot=Hxtot_d if abs(Hxtot_d) > Blimit else 0

    H1y=H1y_d if abs(H1y_d) > Blimit else 0
    H2y=H2y_d if abs(H2y_d) > Blimit else 0
    Hy_d=H1y_d+H2y_d
    Hy=Hy_d if abs(Hy_d) > Blimit else 0
    H1Cy=H1Cy_d if abs(H1Cy_d) > Blimit else 0
    H2Cy=H2Cy_d if abs(H2Cy_d) > Blimit else 0
    HCy_d=H1Cy_d+H2Cy_d
    HCy=HCy_d if abs(HCy_d) > Blimit else 0
    if(Btype==0):
      Hytot_d=H1y_d+H2y_d+H1Cy_d+H2Cy_d
    elif(Btype==1):
      Hytot_d=Hy
    elif(Btype==2):
      Hytot_d=Hyc
      
    Hytot=Hytot_d if abs(Hytot_d) > Blimit else 0
 
    H1Cz=H1Cz_d if abs(H1Cz_d) > Blimit else 0
    H2Cz=H2Cz_d if abs(H2Cz_d) > Blimit else 0
    if(Btype!=1):
      Hz_d=H1Cz_d+H2Cz_d
    else:
      Hz_d=0
    Hztot=Hz_d if abs(Hz_d) > Blimit else 0

    dateoe.write('{:14.8e}'.format(Hxtot)+sp+'{:14.8e}'.format(Hytot)+sp+'{:14.8e}'.format(Hztot)+"\n")
    dateog.write(x+sp+y+sp+z+sp+'{:14.8e}'.format(H1x*convfactor)+sp+'{:14.8e}'.format(H2x*convfactor)+sp+'{:14.8e}'.format(Hx*convfactor)+sp+'{:14.8e}'.format(H1Cx*convfactor)+sp+'{:14.8e}'.format(H2Cx*convfactor)+sp+'{:14.8e}'.format(HCx*convfactor)+sp+'{:14.8e}'.format(Hxtot*convfactor)+sp+'{:14.8e}'.format(H1y*convfactor)+sp+'{:14.8e}'.format(H2y*convfactor)+sp+'{:14.8e}'.format(Hy*convfactor)+sp+'{:14.8e}'.format(H1Cy*convfactor)+sp+'{:14.8e}'.format(H2Cy*convfactor)+sp+'{:14.8e}'.format(HCy*convfactor)+sp+'{:14.8e}'.format(Hytot*convfactor)+sp+'{:14.8e}'.format(H1Cz*convfactor)+sp+'{:14.8e}'.format(H2Cz*convfactor)+sp+'{:14.8e}'.format(Hztot*convfactor)+"\n")
  else: #we just report the blank line in the output file
    dateog.write("\n")
 

dateoe.close()
dateog.close()
print("All done.\n")
