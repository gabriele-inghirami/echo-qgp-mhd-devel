# *****************************************************************************
# *                                                                           *
# *  ECHO-QGP                                                                 *         
# *                                                                           *         
# *  Version: 1.5.0-alpha                                                     *
# *                                                                           *
# *  Copyright (C) 2018-2023 The ECHO-QGP team                                * 
# *                                                                           *
# *  File: dan.py                                                             *
# *                                                                           *
# *  License: GPL version 2.0 (Please, read the file LICENSE.TXT)             *
# *                                                                           *
# *  This program is free software; you can redistribute it and/or            *
# *  modify it under the terms of the GNU General Public License              *
# *  as published by the Free Software Foundation; either version 2           *
# *  of the License, or (at your option) any later version.                   *
# *                                                                           *
# *  This program is distributed in the hope that it will be useful,          *
# *  but WITHOUT ANY WARRANTY; without even the implied warranty of           *
# *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
# *  GNU General Public License for more details.                             *
# *                                                                           *
# *  You should have received a copy of the GNU General Public License        *
# *  along with this program; if not, write to the Free Software              *
# *  Foundation, Inc., 51 Franklin Street, Fifth Floor,                       *
# *  Boston, MA  02110-1301, USA.                                             *
# *                                                                           *
# *  Authors: Gabriele Inghirami (g.inghirami@gsi.de)                         *
# *                                                                           * 
# *****************************************************************************

import fileinput
import math
import numpy as np
import sys
import os
import errno
import time
import scipy.interpolate as intp

var_list = []
t = x = y = z = 0.
rho = vx = vy = vz = en = pr = temp = entr = g = tfo = 0.
bulk = pitt = pitx = pity = pitz = pixy = pixz = piyz = pixx = piyy = pizz = 0.
bx = by = bz = ex = ey = ez = glm = rc = 0.
xmin = xmax = ymin = ymax = zmin = zmax = 0.
xstep = ystep = zstep = 0.
nx = ny = nz = 0
coord = 0
prec = 8
vis = False
mhd = False

def load(infile_name,*args):
  """It loads the data contained in the output files and stores them into numpy arrays.\nSyntax: load("datafile","var1","var2"...)\nInput arguments:\n"datafile" : the name of the ECHO-QGP output file to read, e.g. "data0003.dat"\n
"var1", "var2"...: (optional) The name of the variables to be read (e.g. "vx" or "en"). If no name is given, all available variables are read by default.\nReturn value: 0 in case of success.\nExample: dan.load("data0001.dat","vx","vy","vz")\n"""

  global x,y,z,rho,vx,vy,vz,en,pr,temp,entr,bulk,pitt,pitx,pity,pitz,pixx,pixy,pixz,piyz,piyy,pizz,g,bx,by,bz,ex,ey,ez,glm,rc
  global t,nx,ny,nz,ncells,xmin,xmax,ymin,ymax,coord,prec,vis,mhd,tfo,xstep,ystep,zstep,var_list,bimp,loaded,otype
  
  loaded=False #it will become True at the end of the function, if reached, and returned to the caller
  if(infile_name==''):
     print("Input file name missing.\n")
     print("Syntax: python3 eload.py <inputfile> <list of the variables to real, e.g. vx vy vz. No list is equivalent to include all of them.\n")
     sys.exit(1)
  
  if(len(args)==0):
     print("No list of variables given, I assume that you want to use all the available ones")
     allvars=True 
  else:
     allvars=False 

  #we check that we have the essential file config_summary.dat
  try:
    fp=open("config_summary.dat","r")
  except FileNotFoundError:
    try:
      dnm=os.path.dirname(infile_name)
      fp=open(dnm+"/config_summary.dat","r")
    except FileNotFoundError:
      print("Sorry, I cannot find the file 'config_summary.dat'. (It must be located in the current working directory or in the same directory as config_summary.dat)\n")
      sys.exit(1)

  #we store the informations in a string and we close the file  
  confdata=fp.read()
  fp.close()


  #We determine how large is the grid:
  ss=confdata.find("NX.......=")
  if(ss != -1):
    ee=confdata.find("\n",ss)
    nx=int(confdata[ss+10:ee].strip())
  else:
    print("Error, I cannot determine NX\n")
    sys.exit(1)

  ss=confdata.find("NY.......=")
  if(ss != -1):
    ee=confdata.find("\n",ss)
    ny=int(confdata[ss+10:ee].strip())
  else:
    print("Error, I cannot determine NY\n")
    sys.exit(1)

  ss=confdata.find("NZ.......=")
  if(ss != -1):
    ee=confdata.find("\n",ss)
    nz=int(confdata[ss+10:ee].strip())
  else:
    print("Error, I cannot determine NZ\n")
    sys.exit(1)

  ncells=nx*ny*nz

  ss=confdata.find("XMIN.....=")
  if(ss != -1):
    ee=confdata.find("\n",ss)
    xmin=np.float64(confdata[ss+10:ee].strip())
  else:
    print("Error, I cannot determine XMIN\n")
    sys.exit(1)

  ss=confdata.find("XMAX.....=")
  if(ss != -1):
    ee=confdata.find("\n",ss)
    xmax=np.float64(confdata[ss+10:ee].strip())
  else:
    print("Error, I cannot determine XMAX\n")
    sys.exit(1)
  
  ss=confdata.find("YMIN.....=")
  if(ss != -1):
    ee=confdata.find("\n",ss)
    ymin=np.float64(confdata[ss+10:ee].strip())
  else:
    print("Error, I cannot determine YMIN\n")
    sys.exit(1)

  ss=confdata.find("YMAX.....=")
  if(ss != -1):
    ee=confdata.find("\n",ss)
    ymax=np.float64(confdata[ss+10:ee].strip())
  else:
    print("Error, I cannot determine YMAX\n")
    sys.exit(1)

  ss=confdata.find("ZMIN.....=")
  if(ss != -1):
    ee=confdata.find("\n",ss)
    zmin=np.float64(confdata[ss+10:ee].strip())
  else:
    print("Error, I cannot determine ZMIN\n")
    sys.exit(1)

  ss=confdata.find("ZMAX.....=")
  if(ss != -1):
    ee=confdata.find("\n",ss)
    zmax=np.float64(confdata[ss+10:ee].strip())
  else:
    print("Error, I cannot determine ZMAX\n")
    sys.exit(1)

  #now we get informations about the type of simulation
  ss=confdata.find("COORD....=")
  if(ss != -1):
    ee=confdata.find("\n",ss)
    coord=int(confdata[ss+10:ee].strip())
  else:
    print("Error, I cannot determine COORD, i.e. the coordinate type\n")
    sys.exit(1)
  #we check that the info about the coordinate type makes sens
  if(coord==1):
    "Data about simulation in Minkowski coordinates"
  elif(coord==2):
    "Data about simulation in Milne/Bjorken coordinates"
  else:
    "Data about simulation in unknown coordinate system... I quit."
    sys.exit(2)

  ss=confdata.find("B........=")
  if(ss != -1):
    ee=confdata.find("\n",ss)
    bimpstring=confdata[ss+10:ee].strip()
    bimp=np.float64(bimpstring)
  else:
    print("Error, I cannot determine which was the impact parameter b. I quit.\n")
    sys.exit(1)

  ss=confdata.find("OUTP_PREC=")
  if(ss != -1):
    ee=confdata.find("\n",ss)
    prec=confdata[ss+10:ee].strip()
    if(prec == "8"):
      otype=np.float64
      dimfloat=8
      print("Output precision: 8 bytes (double)\n")
    else:
      otype=np.float32
      dimfloat=4
      print("Output precision: 4 bytes (single)\n")
  else:
    print("Error, I cannot determine the output precision (i.e. 8 or 4 bytes)\n")
    sys.exit(1)

  ss=confdata.find("VISCOUS..=")
  if(ss != -1):
    ee=confdata.find("\n",ss)
    vis_flag=confdata[ss+10:ee].strip()
    if(vis_flag == "1"):
      vis=True
    else:
      vis=False
  else:
    print("Error, I cannot determine VISCOUS, i.e. if the simulation is viscous or not\n")
    sys.exit(1)

  ss=confdata.find("MHD......=")
  if(ss != -1):
    ee=confdata.find("\n",ss)
    mhd_flag=confdata[ss+10:ee].strip()
    if(mhd_flag == "1"):
      mhd=True
    else:
      mhd=False
  else:
    print("Error, I cannot determine MHD\n")
    sys.exit(1)

  ss=confdata.find("FREEZKIND=")
  if(ss != -1):
    ee=confdata.find("\n",ss)
    freezekind=int(confdata[ss+10:ee].strip())
    if(freezekind == 0):
      print("Freeze-out hypersurface based on temperature\n")
      hyp=True
    elif(freezekind == 1):
      print("Freeze-out hypersurface based on energy density\n")
      hyp=True
    else:
      print("I cannot determine how the freeze-out hypersurface is determined, allowed values are only 0 or 1 and I got "+str(freezekind)+"\n")
      sys.exit(1)
  else:
    print("I cannot determine how the freeze-out hypersurface is determined, I assume that it was not built.\n")
    hyp=False

  if(hyp):
    ss=confdata.find("FREEZEVAL=")
    if(ss != -1):
      ee=confdata.find("\n",ss)
      tfo=otype(confdata[ss+10:ee].strip())
    else:
      print("I did not find the value of the freezeout temperature, I assume that it is 0.\n")
      tfo=0.
  else:
    tfo=0.


  if(mhd and vis):
    print("Sorry, but at the moment viscous and MHD modes are incompatible...\nThere is something inconsistent in your config_summary.dat file\n")
    sys.exit(2)

  #now we want to know which variables are available
  ind_cur=dimfloat #index pointer of the output file, we start after the number of bytes of the first entry, which is the time
  blockbytes=ncells*dimfloat


  ss=confdata.find("density..=")
  if(ss != -1):
    ee=confdata.find("\n",ss)
    rho_flag=confdata[ss+10:ee].strip()
    if(rho_flag == "1"):
      rho_printed=True
      ind_rho=ind_cur
      ind_cur=ind_cur+blockbytes
    else:
      rho_printed=False
  else:
    print("Error, I cannot determine whether the density was included in the output or not\n")
    sys.exit(1)

  ss=confdata.find("vx.......=")
  if(ss != -1):
    ee=confdata.find("\n",ss)
    vx_flag=confdata[ss+10:ee].strip()
    if(vx_flag == "1"):
      vx_printed=True
      ind_vx=ind_cur
      ind_cur=ind_cur+blockbytes
    else:
      vx_printed=False
  else:
    print("Error, I cannot determine whether vx was included in the output or not\n")
    sys.exit(1)

  ss=confdata.find("vy.......=")
  if(ss != -1):
    ee=confdata.find("\n",ss)
    vy_flag=confdata[ss+10:ee].strip()
    if(vy_flag == "1"):
      vy_printed=True
      ind_vy=ind_cur
      ind_cur=ind_cur+blockbytes
    else:
      vy_printed=False
  else:
    print("Error, I cannot determine whether vy was included in the output or not\n")
    sys.exit(1)

  ss=confdata.find("vz.......=")
  if(ss != -1):
    ee=confdata.find("\n",ss)
    vz_flag=confdata[ss+10:ee].strip()
    if(vz_flag == "1"):
      vz_printed=True
      ind_vz=ind_cur
      ind_cur=ind_cur+blockbytes
    else:
      vz_printed=False
  else:
    print("Error, I cannot determine whether vz was included in the output or not\n")
    sys.exit(1)

  ss=confdata.find("pressure.=")
  if(ss != -1):
    ee=confdata.find("\n",ss)
    pr_flag=confdata[ss+10:ee].strip()
    if(pr_flag == "1"):
      pr_printed=True
      ind_pr=ind_cur
      ind_cur=ind_cur+blockbytes
    else:
      pr_printed=False
  else:
    print("Error, I cannot determine whether the pressure was included in the output or not\n")
    sys.exit(1)

  ss=confdata.find("ene_dens.=")
  if(ss != -1):
    ee=confdata.find("\n",ss)
    en_flag=confdata[ss+10:ee].strip()
    if(en_flag == "1"):
      en_printed=True
      ind_en=ind_cur
      ind_cur=ind_cur+blockbytes
    else:
      en_printed=False
  else:
    print("Error, I cannot determine whether the energy density was included in the output or not\n")
    sys.exit(1)

  ss=confdata.find("temper...=")
  if(ss != -1):
    ee=confdata.find("\n",ss)
    temp_flag=confdata[ss+10:ee].strip()
    if(temp_flag == "1"):
      temp_printed=True
      ind_temp=ind_cur
      ind_cur=ind_cur+blockbytes
    else:
      temp_printed=False
  else:
    print("Error, I cannot determine whether the temperature was included in the output or not\n")
    sys.exit(1)

  ss=confdata.find("entr_dens=")
  if(ss != -1):
    ee=confdata.find("\n",ss)
    entr_flag=confdata[ss+10:ee].strip()
    if(entr_flag == "1"):
      entr_printed=True
      ind_entr=ind_cur
      ind_cur=ind_cur+blockbytes
    else:
      entr_printed=False
  else:
    print("Error, I cannot determine whether the entropy density was included in the output or not\n")
    sys.exit(1)

  if(vis):
    ss=confdata.find("bulk_visc=")
    if(ss != -1):
      ee=confdata.find("\n",ss)
      bulk_flag=confdata[ss+10:ee].strip()
      if(bulk_flag == "1"):
        bulk_printed=True
        ind_bulk=ind_cur
        ind_cur=ind_cur+blockbytes
      else:
        bulk_printed=False
    else:
      print("Error, I cannot determine whether the bulk viscosity was included in the output or not\n")
      sys.exit(1)

    ss=confdata.find("pi^tt....=")
    if(ss != -1):
      ee=confdata.find("\n",ss)
      pitt_flag=confdata[ss+10:ee].strip()
      if(pitt_flag == "1"):
        pitt_printed=True
        ind_pitt=ind_cur
        ind_cur=ind_cur+blockbytes
      else:
        pitt_printed=False
    else:
      print("Error, I cannot determine whether the pi^tt shear tensor component was included in the output or not\n")
      sys.exit(1)

    ss=confdata.find("pi^tx....=")
    if(ss != -1):
      ee=confdata.find("\n",ss)
      pitx_flag=confdata[ss+10:ee].strip()
      if(pitx_flag == "1"):
        pitx_printed=True
        ind_pitx=ind_cur
        ind_cur=ind_cur+blockbytes
      else:
        pitx_printed=False
    else:
      print("Error, I cannot determine whether the pi^tx shear tensor component was included in the output or not\n")
      sys.exit(1)

    ss=confdata.find("pi^ty....=")
    if(ss != -1):
      ee=confdata.find("\n",ss)
      pity_flag=confdata[ss+10:ee].strip()
      if(pity_flag == "1"):
        pity_printed=True
        ind_pity=ind_cur
        ind_cur=ind_cur+blockbytes
      else:
        pity_printed=False
    else:
      print("Error, I cannot determine whether the pi^ty shear tensor component was included in the output or not\n")
      sys.exit(1)

    ss=confdata.find("pi^tz....=")
    if(ss != -1):
      ee=confdata.find("\n",ss)
      pitz_flag=confdata[ss+10:ee].strip()
      if(pitz_flag == "1"):
        pitz_printed=True
        ind_pitz=ind_cur
        ind_cur=ind_cur+blockbytes
      else:
        pitz_printed=False
    else:
      print("Error, I cannot determine whether the pi^tz shear tensor component was included in the output or not\n")
      sys.exit(1)

    ss=confdata.find("pi^xy....=")
    if(ss != -1):
      ee=confdata.find("\n",ss)
      pixy_flag=confdata[ss+10:ee].strip()
      if(pixy_flag == "1"):
        pixy_printed=True
        ind_pixy=ind_cur
        ind_cur=ind_cur+blockbytes
      else:
        pixy_printed=False
    else:
      print("Error, I cannot determine whether the pi^xy shear tensor component was included in the output or not\n")
      sys.exit(1)

    ss=confdata.find("pi^xz....=")
    if(ss != -1):
      ee=confdata.find("\n",ss)
      pixz_flag=confdata[ss+10:ee].strip()
      if(pixz_flag == "1"):
        pixz_printed=True
        ind_pixz=ind_cur
        ind_cur=ind_cur+blockbytes
      else:
        pixz_printed=False
    else:
      print("Error, I cannot determine whether the pi^xz shear tensor component was included in the output or not\n")
      sys.exit(1)

    ss=confdata.find("pi^yz....=")
    if(ss != -1):
      ee=confdata.find("\n",ss)
      piyz_flag=confdata[ss+10:ee].strip()
      if(piyz_flag == "1"):
        piyz_printed=True
        ind_piyz=ind_cur
        ind_cur=ind_cur+blockbytes
      else:
        piyz_printed=False
    else:
      print("Error, I cannot determine whether the pi^yz shear tensor component was included in the output or not\n")
      sys.exit(1)

    ss=confdata.find("pi^xx....=")
    if(ss != -1):
      ee=confdata.find("\n",ss)
      pixx_flag=confdata[ss+10:ee].strip()
      if(pixx_flag == "1"):
        pixx_printed=True
        ind_pixx=ind_cur
        ind_cur=ind_cur+blockbytes
      else:
        pixx_printed=False
    else:
      print("Error, I cannot determine whether the pi^xx shear tensor component was included in the output or not\n")
      sys.exit(1)

    ss=confdata.find("pi^yy....=")
    if(ss != -1):
      ee=confdata.find("\n",ss)
      piyy_flag=confdata[ss+10:ee].strip()
      if(piyy_flag == "1"):
        piyy_printed=True
        ind_piyy=ind_cur
        ind_cur=ind_cur+blockbytes
      else:
        piyy_printed=False
    else:
      print("Error, I cannot determine whether the pi^yy shear tensor component was included in the output or not\n")
      sys.exit(1)

    ss=confdata.find("pi^zz....=")
    if(ss != -1):
      ee=confdata.find("\n",ss)
      pizz_flag=confdata[ss+10:ee].strip()
      if(pizz_flag == "1"):
        pizz_printed=True
        ind_pizz=ind_cur
        ind_cur=ind_cur+blockbytes
      else:
        pizz_printed=False
    else:
      print("Error, I cannot determine whether the pi^zz shear tensor component was included in the output or not\n")
      sys.exit(1)

  ss=confdata.find("gamma....=")
  if(ss != -1):
    ee=confdata.find("\n",ss)
    gamma_flag=confdata[ss+10:ee].strip()
    if(gamma_flag == "1"):
      gamma_printed=True
      ind_gamma=ind_cur
      ind_cur=ind_cur+blockbytes
    else:
      gamma_printed=False
  else:
    print("Error, I cannot determine whether the gamma Lorentz factor was included in the output or not\n")
    sys.exit(1)

  if(mhd):
    ss=confdata.find("bx.......=")
    if(ss != -1):
      ee=confdata.find("\n",ss)
      bx_flag=confdata[ss+10:ee].strip()
      if(bx_flag == "1"):
        bx_printed=True
        ind_bx=ind_cur
        ind_cur=ind_cur+blockbytes
      else:
        bx_printed=False
    else:
      print("Error, I cannot determine whether bx was included in the output or not\n")
      sys.exit(1)

    ss=confdata.find("by.......=")
    if(ss != -1):
      ee=confdata.find("\n",ss)
      by_flag=confdata[ss+10:ee].strip()
      if(by_flag == "1"):
        by_printed=True
        ind_by=ind_cur
        ind_cur=ind_cur+blockbytes
      else:
        by_printed=False
    else:
      print("Error, I cannot determine whether by was included in the output or not\n")
      sys.exit(1)

    ss=confdata.find("bz.......=")
    if(ss != -1):
      ee=confdata.find("\n",ss)
      bz_flag=confdata[ss+10:ee].strip()
      if(bz_flag == "1"):
        bz_printed=True
        ind_bz=ind_cur
        ind_cur=ind_cur+blockbytes
      else:
        bz_printed=False
    else:
      print("Error, I cannot determine whether bz was included in the output or not\n")
      sys.exit(1)

    ss=confdata.find("ex.......=")
    if(ss != -1):
      ee=confdata.find("\n",ss)
      ex_flag=confdata[ss+10:ee].strip()
      if(ex_flag == "1"):
        ex_printed=True
        ind_ex=ind_cur
        ind_cur=ind_cur+blockbytes
      else:
        ex_printed=False
    else:
      print("Error, I cannot determine whether ex was included in the output or not\n")
      sys.exit(1)

    ss=confdata.find("ey.......=")
    if(ss != -1):
      ee=confdata.find("\n",ss)
      ey_flag=confdata[ss+10:ee].strip()
      if(ey_flag == "1"):
        ey_printed=True
        ind_ey=ind_cur
        ind_cur=ind_cur+blockbytes
      else:
        ey_printed=False
    else:
      print("Error, I cannot determine whether ey was included in the output or not\n")
      sys.exit(1)

    ss=confdata.find("ez.......=")
    if(ss != -1):
      ee=confdata.find("\n",ss)
      ez_flag=confdata[ss+10:ee].strip()
      if(ez_flag == "1"):
        ez_printed=True
        ind_ez=ind_cur
        ind_cur=ind_cur+blockbytes
      else:
        ez_printed=False
    else:
      print("Error, I cannot determine whether ez was included in the output or not\n")
      sys.exit(1)

    ss=confdata.find("glm......=")
    if(ss != -1):
      ee=confdata.find("\n",ss)
      glm_flag=confdata[ss+10:ee].strip()
      if(glm_flag == "1"):
        glm_printed=True
        ind_glm=ind_cur
        ind_cur=ind_cur+blockbytes
      else:
        glm_printed=False
    else:
      print("Error, I cannot determine whether glm was included in the output or not\n")
      sys.exit(1)

    ss=confdata.find("rc.......=")
    if(ss != -1):
      ee=confdata.find("\n",ss)
      rc_flag=confdata[ss+10:ee].strip()
      if(rc_flag == "1"):
        rc_printed=True
        ind_rc=ind_cur
        ind_cur=ind_cur+blockbytes
      else:
        rc_printed=False
    else:
      print("Error, I cannot determine whether rc (the electric charge in the comoving frame) was included in the output or not\n")
      sys.exit(1)


  #now we use the informations collected so far
  #first, we build the grid
  xstep=(xmax-xmin)/nx
  ystep=(ymax-ymin)/ny
  zstep=(zmax-zmin)/nz
  x=np.linspace(xmin+xstep/2,xmax-xstep/2,nx,endpoint=True,dtype=np.float64)
  y=np.linspace(ymin+ystep/2,ymax-ystep/2,ny,endpoint=True,dtype=np.float64)
  z=np.linspace(zmin+zstep/2,zmax-zstep/2,nz,endpoint=True,dtype=np.float64)

  inline=' '.join(args[:])

  try:
    infile=open(infile_name,"r+b")
  except FileNotFoundError:
    print(infile_name+" cannot be opened. I quit.\n")
    sys.exit(2)
  
  var_list=[]

  t=np.fromfile(infile,dtype=otype,count=1)[0]
  print("Time is stored in the dan.t float variable\n") 
  
  if(allvars or (inline.find("rho") != -1)):
    if(rho_printed):
      infile.seek(ind_rho, os.SEEK_SET)
      rho=np.fromfile(infile,dtype=otype,count=ncells).reshape([nz,ny,nx]).transpose()
      print("Density stored in the dan.rho array\n") 
      var_list.append('rho')
    else:
      print("Density not recorded\n") 

  if(allvars or (inline.find("vx") != -1)):
    if(vx_printed):
      infile.seek(ind_vx, os.SEEK_SET)
      vx=np.fromfile(infile,dtype=otype,count=ncells).reshape([nz,ny,nx]).transpose() 
      print("vx stored in the dan.vx array\n") 
      var_list.append('vx')
    else:
      print("vx not recorded\n") 

  if(allvars or (inline.find("vy") != -1)):
    if(vx_printed):
      infile.seek(ind_vy, os.SEEK_SET)
      vy=np.fromfile(infile,dtype=otype,count=ncells).reshape([nz,ny,nx]).transpose() 
      print("vy stored in the dan.vy array\n") 
      var_list.append('vy')
    else:
      print("vy not recorded\n") 

  if(allvars or (inline.find("vz") != -1)):
    if(vz_printed):
      infile.seek(ind_vz, os.SEEK_SET)
      vz=np.fromfile(infile,dtype=otype,count=ncells).reshape([nz,ny,nx]).transpose() 
      print("vz stored in the dan.vz array\n") 
      var_list.append('vz')
    else:
      print("vz not recorded\n") 

  if(allvars or (inline.find("pr") != -1)):
    if(pr_printed):
      infile.seek(ind_pr, os.SEEK_SET)
      pr=np.fromfile(infile,dtype=otype,count=ncells).reshape([nz,ny,nx]).transpose() 
      print("pressure stored in the dan.pr array\n") 
      var_list.append('pr')
    else:
      print("pressure not recorded\n") 

  if(allvars or (inline.find("en") != -1)):
    if(en_printed):
      infile.seek(ind_en, os.SEEK_SET)
      en=np.fromfile(infile,dtype=otype,count=ncells).reshape([nz,ny,nx]).transpose() 
      print("energy density stored in the dan.en array\n") 
      var_list.append('en')
    else:
      print("energy density not recorded\n") 
  
  if(allvars or (inline.find("temp") != -1)):
    if(temp_printed):
      infile.seek(ind_temp, os.SEEK_SET)
      temp=np.fromfile(infile,dtype=otype,count=ncells).reshape([nz,ny,nx]).transpose() 
      print("temperature stored in the dan.temp array\n") 
      var_list.append('temp')
    else:
      print("temperature not recorded\n") 
  
  if(allvars or (inline.find("entr") != -1)):
    if(entr_printed):
      infile.seek(ind_entr, os.SEEK_SET)
      entr=np.fromfile(infile,dtype=otype,count=ncells).reshape([nz,ny,nx]).transpose() 
      print("entropy density stored in the dan.entr array\n") 
      var_list.append('entr')
    else:
      print("entropy density not recorded\n") 
  
  if(vis):
  
    if(allvars or (inline.find("bulk") != -1)):
      if(bulk_printed):
        infile.seek(ind_bulk, os.SEEK_SET)
        bulk=np.fromfile(infile,dtype=otype,count=ncells).reshape([nz,ny,nx]).transpose() 
        print("bulk viscosity stored in the dan.bulk array\n") 
        var_list.append('bulk')
      else:
        print("bulk viscosity not recorded\n") 
  
    if(allvars or (inline.find("pitt") != -1)):
      if(pitt_printed):
        infile.seek(ind_pitt, os.SEEK_SET)
        pitt=np.fromfile(infile,dtype=otype,count=ncells).reshape([nz,ny,nx]).transpose() 
        print("pi^tt stored in the dan.pitt array\n") 
        var_list.append('pitt')
      else:
        print("pitt not recorded\n") 
  
    if(allvars or (inline.find("pitx") != -1)):
      if(pitx_printed):
        infile.seek(ind_pitx, os.SEEK_SET)
        pitx=np.fromfile(infile,dtype=otype,count=ncells).reshape([nz,ny,nx]).transpose() 
        print("pi^tx stored in the dan.pitx array\n") 
        var_list.append('pitx')
      else:
        print("pitx not recorded\n") 
  
    if(allvars or (inline.find("pity") != -1)):
      if(pity_printed):
        infile.seek(ind_pity, os.SEEK_SET)
        pity=np.fromfile(infile,dtype=otype,count=ncells).reshape([nz,ny,nx]).transpose() 
        print("pi^ty stored in the dan.pity array\n") 
        var_list.append('pity')
      else:
        print("pity not recorded\n") 
  
    if(allvars or (inline.find("pitz") != -1)):
      if(pitz_printed):
        infile.seek(ind_pitz, os.SEEK_SET)
        pitz=np.fromfile(infile,dtype=otype,count=ncells).reshape([nz,ny,nx]).transpose() 
        print("pi^tz stored in the dan.pitz array\n") 
        var_list.append('pitz')
      else:
        print("pitz not recorded\n") 
  
    if(allvars or (inline.find("pixy") != -1)):
      if(pixy_printed):
        infile.seek(ind_pixy, os.SEEK_SET)
        pixy=np.fromfile(infile,dtype=otype,count=ncells).reshape([nz,ny,nx]).transpose() 
        print("pi^xy stored in the dan.pixy array\n") 
        var_list.append('pixy')
      else:
        print("pixy not recorded\n") 
  
    if(allvars or (inline.find("pixz") != -1)):
      if(pixz_printed):
        infile.seek(ind_pixz, os.SEEK_SET)
        pixz=np.fromfile(infile,dtype=otype,count=ncells).reshape([nz,ny,nx]).transpose() 
        print("pi^xz stored in the dan.pixz array\n") 
        var_list.append('pixz')
      else:
        print("pixz not recorded\n") 
  
    if(allvars or (inline.find("piyz") != -1)):
      if(piyz_printed):
        infile.seek(ind_piyz, os.SEEK_SET)
        piyz=np.fromfile(infile,dtype=otype,count=ncells).reshape([nz,ny,nx]).transpose() 
        print("pi^yz stored in the dan.piyz array\n") 
        var_list.append('piyz')
      else:
        print("piyz not recorded\n") 
  
    if(allvars or (inline.find("pixx") != -1)):
      if(pixx_printed):
        infile.seek(ind_pixx, os.SEEK_SET)
        pixx=np.fromfile(infile,dtype=otype,count=ncells).reshape([nz,ny,nx]).transpose() 
        print("pi^xx stored in the dan.pixx array\n") 
        var_list.append('pixx')
      else:
        print("pixx not recorded\n") 
  
    if(allvars or (inline.find("piyy") != -1)):
      if(piyy_printed):
        infile.seek(ind_piyy, os.SEEK_SET)
        piyy=np.fromfile(infile,dtype=otype,count=ncells).reshape([nz,ny,nx]).transpose() 
        print("pi^yy stored in the dan.piyy array\n") 
        var_list.append('piyy')
      else:
        print("piyy not recorded\n") 
  
    if(allvars or (inline.find("pizz") != -1)):
      if(pizz_printed):
        infile.seek(ind_pizz, os.SEEK_SET)
        pizz=np.fromfile(infile,dtype=otype,count=ncells).reshape([nz,ny,nx]).transpose() 
        print("pi^zz stored in the dan.pizz array\n") 
        var_list.append('pizz')
      else:
        print("pizz not recorded\n") 
  
  if(allvars or (inline.find("gamma") != -1)):
    if(gamma_printed):
      infile.seek(ind_gamma, os.SEEK_SET)
      g=np.fromfile(infile,dtype=otype,count=ncells).reshape([nz,ny,nx]).transpose() 
      print("gamma Lorentz factor stored in the dan.g array\n") 
      var_list.append('gamma')
    else:
      print("gamma Lorentz factor not recorded\n") 
  
  
  if(mhd):
  
    if(allvars or (inline.find("bx") != -1)):
      if(bx_printed):
        infile.seek(ind_bx, os.SEEK_SET)
        bx=np.fromfile(infile,dtype=otype,count=ncells).reshape([nz,ny,nx]).transpose() 
        print("bx stored in the dan.bx array\n") 
        var_list.append('bx')
      else:
        print("bx not recorded\n") 
  
    if(allvars or (inline.find("by") != -1)):
      if(by_printed):
        infile.seek(ind_by, os.SEEK_SET)
        by=np.fromfile(infile,dtype=otype,count=ncells).reshape([nz,ny,nx]).transpose() 
        print("by stored in the dan.by array\n") 
        var_list.append('by')
      else:
        print("by not recorded\n") 
  
    if(allvars or (inline.find("bz") != -1)):
      if(bz_printed):
        infile.seek(ind_bz, os.SEEK_SET)
        bz=np.fromfile(infile,dtype=otype,count=ncells).reshape([nz,ny,nx]).transpose() 
        print("bz stored in the dan.bz array\n") 
        var_list.append('bz')
      else:
        print("bz not recorded\n") 
  
    if(allvars or (inline.find("ex") != -1)):
      if(ex_printed):
        infile.seek(ind_ex, os.SEEK_SET)
        ex=np.fromfile(infile,dtype=otype,count=ncells).reshape([nz,ny,nx]).transpose() 
        print("ex stored in the dan.ex array\n") 
        var_list.append('ex')
      else:
        print("ex not recorded\n") 
  
    if(allvars or (inline.find("ey") != -1)):
      if(ey_printed):
        infile.seek(ind_ey, os.SEEK_SET)
        ey=np.fromfile(infile,dtype=otype,count=ncells).reshape([nz,ny,nx]).transpose() 
        print("ey stored in the dan.ey array\n") 
        var_list.append('ey')
      else:
        print("ey not recorded\n") 
  
    if(allvars or (inline.find("ez") != -1)):
      if(ez_printed):
        infile.seek(ind_ez, os.SEEK_SET)
        ez=np.fromfile(infile,dtype=otype,count=ncells).reshape([nz,ny,nx]).transpose() 
        print("ez stored in the dan.ez array\n") 
        var_list.append('ez')
      else:
        print("ez not recorded\n") 
  
    if(allvars or (inline.find("glm") != -1)):
      if(glm_printed):
        infile.seek(ind_glm, os.SEEK_SET)
        glm=np.fromfile(infile,dtype=otype,count=ncells).reshape([nz,ny,nx]).transpose() 
        print("glm stored in the dan.glm array\n") 
        var_list.append('glm')
      else:
        print("glm not recorded\n") 

    if(allvars or (inline.find("rc") != -1)):
      if(rc_printed):
        infile.seek(ind_rc, os.SEEK_SET)
        rc=np.fromfile(infile,dtype=otype,count=ncells).reshape([nz,ny,nx]).transpose() 
        print("rc stored in the dan.rc array\n") 
        var_list.append('rc')
      else:
        print("rc not recorded\n") 
  
  
  infile.close()
  loaded=True
  print("var_list is a list containing the names of all variable arrays.\n")
  return loaded
  
  

def load_der(infile_name):
  """It loads the data about partial derivatives contained in the output files and stores them into numpy arrays.\nSyntax: load_der("datafile")\nInput arguments:\n"datafile" : the name of the ECHO-QGP output file with partial derivatives to read, e.g. "der0003.dat"\n
Return value: 0 in case of success.\nExample: dan.load_der("der0001.dat")\n"""

  global derivs, dtx, dty, dtz, dxx, dxy, dxz, dyx, dyy, dyz, dzx, dzy, dzz, dtt, dxt, dyt, dzt, dtet, loaded_der
  dtx=0
  dty=1
  dty=2
  dxx=3
  dxy=4
  dxz=5
  dyz=6
  dyy=7
  dyz=8
  dzx=9
  dzy=10
  dzz=11
  dtt=12
  dxt=13
  dyt=14
  dzt=15
  dtet=16
  dtex=17
  dtey=18
  dtez=19
  
  loaded_der=False
  
  if(not loaded):
	  print("Sorry, but it seems that you did not load the file with main variables referring to the same simulation timestep.\nPlease, load it before calling this function")
	  return loaded_der
  
  if(infile_name==''):
     print("Input file name missing.\n")
     print("Syntax: python3 eload.py <inputfile> <list of the variables to real, e.g. vx vy vz. No list is equivalent to include all of them.\n")
     sys.exit(1)
     
  try:
    infile=open(infile_name,"r+b")
  except FileNotFoundError:
    print(infile_name+" cannot be opened. I quit.\n")
    sys.exit(2)
    
  t_der=np.fromfile(infile,dtype=otype,count=1)[0]
  if(t_der!=t):
	  print("Sorry, but it seems that you did not load the file with main variables referring to the same simulation timestep.\nPlease, load it before calling this function")
	  return loaded_der
	  
  derivs=np.fromfile(infile,dtype=otype,count=nx*ny*nz*20).reshape([nz,ny,nx,20]).transpose()
  print("the partial derivatives are stored in the dan.derivs array\n") 
  var_list.append('derivs')
  loaded_der=True
  return loaded_der
  

#useful functions

#this function computes the total energy over the grid in the pure hd case (with Milne coordinates)
def comp_en_hd(tlim=tfo,ped=False):
    """It returns the total energy contained in the grid, expressed in GeV, without taking into account any EM field and assuming Milne coordiantes.\nSyntax: comp_en_hd({tlim},{ped})\nInput arguments:\ntlim : the minimum value of the temperature in GeV to include a cell in the computation (type: float, default value: the freezeout temperature)\nped : wheter to subtract (True) or not (False, default) the minimum value of the energy density (type: logical, default value: False).\nReturn value: the total energy in GeV, type: float.\nExample: comp_en_hd(0.,False)\n"""
    global t, nx, ny, nz, en, pr, g, grid

    E=0.

    if(ped):
      limit_en=np.amin(en)
      limit_pr=np.amin(pr)
    else:
      limit_en=0.
      limit_pr=0.

    for cx in range(0,nx):
      for cy in range(0,ny):
        for cz in range(0,nz):
          if(temp[cx,cy,cz] > tlim):
            en_dens=en[cx,cy,cz]-limit_en
            p=pr[cx,cy,cz]-limit_pr
            h=en_dens+p
            vveta=vz[cx,cy,cz]
            glf=g[cx,cy,cz]
            ueta=glf*vveta
            eta=z[cz]
            u0=math.cosh(eta)*glf+t*math.sinh(eta)*ueta
            uz=math.sinh(eta)*glf+t*math.cosh(eta)*ueta
            t00=h*u0*u0-p
            t0z=h*u0*uz

            E=E+(xstep*ystep*zstep)*t*(math.cosh(eta)*t00-math.sinh(eta)*t0z)
    
    return E
     
#these functions return the closest array index to the given coordinate value
def xx(point=0.):
    """It returns the closest x array index corresponding to the given x value.\nSyntax: xx(<xpoint>)\nInput arguments:\nxpoint : the value of an x coordinate (type: float, default value: 0.)\nReturn value: the index of the x array element which is closest to xpoint (in case of equal distance of two elements, it returns the smallest one), type: integer.\nExample: xx(0.8)\n"""
    return np.abs(x-point).argmin()

def yy(point=0.):
    """It returns the closest y array index corresponding to the given y value.\nSyntax: yy(<ypoint>)\nInput arguments:\nypoint : the value of an y coordinate (type: float, default value: 0.)\nReturn value: the index of the y array element which is closest to ypoint (in case of equal distance of two elements, it returns the smallest one), type: integer.\nExample: yy(0.8)\n"""
    return np.abs(y-point).argmin()

def zz(point=0.):
    """It returns the closest z array index corresponding to the given z value.\nSyntax: zz(<zpoint>)\nInput arguments:\nzpoint : the value of a z coordinate (type: float, default value: 0.)\nReturn value: the index of the z array element which is closest to zpoint (in case of equal distance of two elements, it returns the smallest one), type: integer.\nExample: zz(0.8)\n"""
    return np.abs(z-point).argmin()

#utilities to create output ascii files to be plotted with gnuplot
def print_to_file(datas,outputfile,comment=""):
    """It prints an array to an ascii file suitable to be plotted with gnuplot (or other tools).\nSyntax: print_to_file(<datas>,<outputfile>,{comment})\nInput arguments:\ndatas: a list or a tuple with one 1,2 or 3D array and 1,2 or 3 1D arrays with the coordinates (type: list or tuple of float arrays),\noutputfile: the name of outputfile (type: string),\ncomment: an optional comment that will be inserted at the beginning of the output file\nReturn value: 0 in case of success (and an ascii output file named outputfile), 1 in case of error in the argument list.\nExample: print_to_file((en,x,y,zz(0)),"en_vs_xy_at_z_eq_0.dat")\n"""
    
    sp="      "    

    nargs=len(datas)
    if(nargs < 2):
      print("Sorry, but you need to provide (as a list or a tuple) both the array with the variable of interest and the 1D arrays containing coordinates. I will not print any file.")
      return 1

    inarr=datas[0]
    ndims=len(inarr.shape)
    if(nargs-1 != ndims):
      print("Sorry, but you did not provide the correct number of 1D arrays containing coordinates. I will not print any file.")
      return 1
    
    x=datas[1]
    xa=len(x)
    if(xa != inarr.shape[0]):
      print("Sorry, but there is a mismatching between the number of elements of the first array with coordinates and the first dimension of the variable array. I will not print any file.")
      return 1
    else:
      dim1=True
      dim2=False
      dim3=False

    ya=1
    if(ndims > 1):
      y=datas[2]
      ya=len(y)
      if(ya != inarr.shape[1]):
        print("Sorry, but there is a mismatching between the number of elements of the second array with coordinates and the second dimension of the variable array. I will not print any file.")
        return 1
      else:
        dim1=False
        dim2=True
        dim3=False

    za=1
    if(ndims > 2):
      z=datas[3]
      za=len(z)
      if(za != inarr.shape[2]):
        print("Sorry, but there is a mismatching between the number of elements of the third array with coordinates and the third dimension of the variable array. I will not print any file.")
        return 1
      else:
        dim1=False
        dim2=False
        dim3=True



    fout=open(outputfile,"w+")
    fout.write("#"+comment+"\n")
    for ix in range(0,xa):  
      xp='{:9.6f}'.format(x[ix])+sp
      if(dim1):
        fout.write(xp+'{:16.14e}'.format(inarr[ix])+"\n")
      else:
        for iy in range(0,ya):
          yp='{:9.6f}'.format(y[iy])+sp
          if(dim2):
              fout.write(xp+yp+'{:16.14e}'.format(inarr[ix,iy])+"\n")
          else:
            for iz in range(0,za):  
              zp='{:9.6f}'.format(z[iz])+sp
              fout.write(xp+yp+zp+'{:16.14e}'.format(inarr[ix,iy,iz])+"\n")

          if(za>1):
            fout.write("\n")
      if(ya>1):
        fout.write("\n")

    return 0 

#utilities to create output ascii files for HF propagation
def print_for_HF(index_start,index_end,eta_cut,flow_outputfile,em_outputfile):
    """It prints two ascii files, one with flow data and the other one with EM field data, to be used for HF propagation.\nWe assume to work with standard output files named out0001.dat, out0002.dat, out0003.dat...\nSyntax: print_for_HF(index_start,index_end,eta_cut,flow_outputfile,em_outputfile)\nInput arguments:\n index_start: the index of the output files from which to start\n index_end: the index of the output files at which to end\neta_cut: values for |eta|>eta_cut (or |z|>eta_cut, depending on the coordinates in use) are not included.\nflow_outputfile: the name of outputfile for the bulk flow data type: string,\nem_outputfile: the name of outputfile for the bulk flow data type: string.\nReturn value: 0 in case of success (and an ascii output file named outputfile).\nExample: print_for_HF((1,110,2.5,"flow_dat_for_charms.dat","em_data_for_charms.dat")\nColumns in the flow_outputfile:  time (fm/c), x(fm), y(fm), eta , Temperature (GeV), density (fm^-3), vx, vy, vz.\nColumns in the em_outputfile:  time (fm/c), x(fm), y(fm), eta ,  Ex , Ey, Ez, Bx, By, Bz. The electromagnetic field is expressed in GeV/fm.\n"""
    
    sp="      "   
    conv=math.sqrt(0.197326)   
    fout_flow=open(flow_outputfile,"w+")
    fout_em=open(em_outputfile,"w+")
    ec=float(eta_cut)
    for it in range(int(index_start),int(index_end+1)):
      infile_name="out"+"{:04d}".format(it)+".dat"
      try:
        infile=open(infile_name,"rb")
      except FileNotFoundError:
        print(infile_name+" cannot be opened. I skip it.\n")
        continue
      infile.close()
      print("Reading "+infile_name+"\n")
      dan.load(infile_name,"vx","vy","vz","en","temp","bx","by","bz","ex","ey","ez")
      tp='{:9.6f}'.format(it)+sp
      for ix in range(0,nx):  
        xp='{:9.6f}'.format(x[ix])+sp
        for iy in range(0,ny):
          yp='{:9.6f}'.format(y[iy])+sp
          for iz in range(0,nz):
             if(abs(z[iz])<=ec):    
                zp='{:9.6f}'.format(z[iz])+sp
                Tv='{:14.8e}'.format(temp[ix,iy,iz])+sp  
                env='{:14.8e}'.format(en[ix,iy,iz])+sp  
                vxv='{:14.8e}'.format(vx[ix,iy,iz])+sp  
                vyv='{:14.8e}'.format(vy[ix,iy,iz])+sp  
                vzv='{:14.8e}'.format(vz[ix,iy,iz])  
                bxv='{:14.8e}'.format(bx[ix,iy,iz]*conv)+sp  
                byv='{:14.8e}'.format(by[ix,iy,iz]*conv)+sp  
                bzv='{:14.8e}'.format(bz[ix,iy,iz]*conv)  
                exv='{:14.8e}'.format(ex[ix,iy,iz]*conv)+sp  
                eyv='{:14.8e}'.format(ey[ix,iy,iz]*conv)+sp  
                ezv='{:14.8e}'.format(ez[ix,iy,iz]*conv)+sp  
                fout_flow.write(tp+xp+yp+zp+Tv+env+vxv+vyv+vzv+"\n")
                fout_em.write(tp+xp+yp+zp+exv+eyv+ezv+bxv+byv+bzv+"\n")

    fout_flow.close()
    fout_em.close()
    print(flow_outputfile+" and "+em_outputfile+" written.")
    return 0 

    
def find_arr(key):
    """It simply returns the array corresponding to the given string name."""
    global rho, vx, vy, vz, en, pr, temp, entr, bulk, pitt, pitx, pity, pitz, pixy, pixz, piyz, pixx, piyy, pizz, g, bx, by, bz, ex, ey, ez, glm, rc
 
    if(key=="rho"):
      return rho
    elif(key=="vx"):
      return vx
    elif(key=="vy"):
      return vy
    elif(key=="vz"):
      return vz
    elif(key=="pr"):
      return pr
    elif(key=="en"):
      return en
    elif(key=="temp"):
      return temp
    elif(key=="entr"):
      return entr
    elif(key=="bulk"):
      return bulk
    elif(key=="pitt"):
      return pitt
    elif(key=="pitx"):
      return pitx
    elif(key=="pity"):
      return pity
    elif(key=="pitz"):
      return pitz
    elif(key=="pixy"):
      return pixy
    elif(key=="piyz"):
      return piyz
    elif(key=="pixz"):
      return pixz
    elif(key=="pixx"):
      return pixx
    elif(key=="piyy"):
      return piyy
    elif(key=="pizz"):
      return pizz
    elif(key=="gamma"):
      return g
    elif(key=="bx"):
      return bx
    elif(key=="by"):
      return by
    elif(key=="bz"):
      return bz
    elif(key=="ex"):
      return ex
    elif(key=="ey"):
      return ey
    elif(key=="ez"):
      return ez
    elif(key=="glm"):
      return glm
    elif(key=="rc"):
      return rc
    else:
      print("Array not found.\n")
      return None

#utility to print some stats
def stats():
    """It provides basic statistical informations, i.e. minimum, mean and maximum values of all variable arrays.\nSyntax: stats()\n"""
    global t, var_list

    print("time: "+str(t))
    for key in var_list:
        myarr=find_arr(key)
        print(key+": min: "+str(np.amin(myarr))+", mean: "+str(np.mean(myarr))+", max: "+str(np.amax(myarr)))

    print("\n")

#it creates a set of standard plots
def p3d(pvars=var_list):
    """It creates a set of standard plots for 3D simulations. For each variable are printed the 1D profiles vs the three axes and the contour and surface plots of the three 2D slices, in each case keeping the value of the not plotted coordinate(s) as close as possible to 0. The plots are saved as png files stored into the check_plots subdirectory.\nSyntax: p3d({pvars})\nInput arguments:\npvars: a list or a tuple of strings with the names of the arryas to be plotted (type: list or tuple of strings, default value: var_list i.e. the list with the names of all available arrays)\nReturn value: 0 in case of success, 1 in case of error in creating the output directory, 2 if an array with a number of dimensions different from 3 is passed to the function.\nExample: p3d(("en","pr"))\n"""
    import matplotlib
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm

    cmap_div=cm.seismic
    cmap_pos=cm.plasma

    global t, var_avail, x, y, z

    print("Using matplotlib version: "+matplotlib.__version__)


    ck="check_plots"

    try:
      os.makedirs(ck)
    except OSError as e:
      if e.errno != errno.EEXIST:
         print("Error in creating the "+ck+" directory.\n")
         raise
         return 1

    tstring='{:05.3f}'.format(t)    
    xc_string='{:05.3f}'.format(x[xx(0)])
    yc_string='{:05.3f}'.format(y[yy(0)])
    zc_string='{:05.3f}'.format(z[zz(0)])

    if(isinstance(pvars, str)):
        tmp_pvars={pvars}
        pvars=tuple(tmp_pvars)
        tmp_pvars=None
 
    for key in pvars:
        print("key is: "+key)
        myarr=find_arr(key) 
        if(key=="rho"):
          ylab="r'$\rho'"
          clm=cmap_pos
        elif(key=="vx"):
          ylab="$v^x$"
          clm=cmap_div
        elif(key=="vy"):
          ylab="$v^y$"
          clm=cmap_div
        elif(key=="vz"):
          ylab="$v^z$"
          clm=cmap_div
        elif(key=="pr"):
          ylab="Pressure [$GeV/fm^3$]"
          clm=cmap_pos
        elif(key=="en"):
          ylab="Energy density [$GeV/fm^3$]"
          clm=cmap_pos
        elif(key=="temp"):
          ylab="Temperature [$GeV$]"
          clm=cmap_pos
        elif(key=="entr"):
          ylab="Entropy density [$fm^{-3}$]"
          clm=cmap_pos
        elif(key=="bulk"):
          ylab="r'$\Pi$'"
          clm=cmap_div
        elif(key=="pitt"):
          ylab="r'$\pi^{tt}$'"
          clm=cmap_div
        elif(key=="pitx"):
          ylab="r'$\pi^{tx}$'"
          clm=cmap_div
        elif(key=="pity"):
          ylab="r'$\pi^{ty}$'"
          clm=cmap_div
        elif(key=="pitz"):
          ylab="r'$\pi^{tz}$'"
          clm=cmap_div
        elif(key=="pixy"):
          ylab="r'$\pi^{xy}$'"
          clm=cmap_div
        elif(key=="piyz"):
          ylab="r'$\pi^{yz}$'"
          clm=cmap_div
        elif(key=="pixz"):
          ylab="r'$\pi^{xz}$'"
          clm=cmap_div
        elif(key=="pixx"):
          ylab="r'$\pi^{xx}$'"
          clm=cmap_div
        elif(key=="piyy"):
          ylab="r'$\pi^{yy}$'"
          clm=cmap_div
        elif(key=="pizz"):
          ylab="r'$\pi^{zz}$'"
          clm=cmap_div
        elif(key=="bx"):
          ylab="'$B^x [GeV^{1/2}/fm^{3/2}]$'"
          clm=cmap_div
        elif(key=="by"):
          ylab="'$B^y [GeV^{1/2}/fm^{3/2}]$'"
          clm=cmap_div
        elif(key=="bz"):
          ylab="'$B^z [GeV^{1/2}/fm^{3/2}]$'"
          clm=cmap_div
        elif(key=="ex"):
          ylab="'$E^x [GeV^{1/2}/fm^{3/2}]$'"
          clm=cmap_div
        elif(key=="ey"):
          ylab="'$E^y [GeV^{1/2}/fm^{3/2}]$'"
          clm=cmap_div
        elif(key=="ez"):
          ylab="'$E^z [GeV^{1/2}/fm^{3/2}]$'"
          clm=cmap_div
        elif(key=="glm"):
          ylab="r'$\psi$'"
          clm=cmap_pos
        elif(key=="rc"):
          ylab="$\\rho_{el}$"
          clm=cmap_div
        else:
          plt.ylabel("  ")   

        if(len(myarr.shape)>3): #we evaluate the number of active directions
          print("Hey, your array has more than 3 dimensions... Something went wrong... I quit.")
          return 2
        elif((nx>1) and (ny>1) and (nz>1)): 
            xyf=True
            xzf=True
            yzf=True
            xpf=True
            ypf=True
            zpf=True
        elif(((nx>1) and (ny>1)) or ((nx>1) and (nz>1)) or ((nz>1) and (ny>1))):
          if(nx==1): #x is the not used axis
            xyf=False
            xzf=False
            yzf=True
            xpf=False
            ypf=True
            zpf=True
          elif(ny==1): #y is the not used axis
            xyf=False
            xzf=True
            yzf=False
            xpf=True
            ypf=False
            zpf=True
          else:
            xyf=True
            xzf=False
            yzf=False
            xpf=True
            ypf=True
            zpf=False
        else:
          if(nx>1): #x is the used axis
            xyf=False
            xzf=False
            yzf=False
            xpf=True
            ypf=False
            zpf=False
          elif(ny>1): #y is the used axis
            xyf=False
            xzf=False
            yzf=False
            xpf=False
            ypf=True
            zpf=False
          else:
            xyf=False
            xzf=False
            yzf=False
            xpf=False
            ypf=False
            zpf=True
         
          
        if(xyf):
          fig = plt.figure()
          fig.set_size_inches(7, 7)
          ax = fig.gca(projection = '3d')
          ax.set_xlabel('y [fm]')
          ax.set_ylabel('x [fm]')
          ax.set_zlabel(ylab)
          X2, Y2 = np.meshgrid(y, x)
          plt.title(key+", t="+tstring+", z="+zc_string)
          surfplot=ax.plot_surface(X2, Y2, myarr[:,:,zz(0)], cmap=clm)
          cbar=plt.colorbar(surfplot,shrink=0.75)
#       plt.savefig(ck+"/"+key+"_3D_xy_t_"+tstring+"_old.png",bbox_inches='tight')
          plt.savefig(ck+"/"+key+"_3D_xy_t_"+tstring+".png")
          fig.clf() 
          plt.close()

        if(yzf):
          fig = plt.figure()
          fig.set_size_inches(7, 7)
          ax = fig.gca(projection = '3d')
          if(coord==1):
            ax.set_xlabel('z')
          else:
            ax.set_xlabel(r'$\eta$')
          ax.set_ylabel('y [fm]')
          ax.set_zlabel(ylab)
          Y2, Z2 = np.meshgrid(z, y)
          plt.title(key+", t="+tstring+", x="+xc_string)
          surfplot=ax.plot_surface(Y2, Z2, myarr[xx(0),:,:], cmap=clm)
          cbar=plt.colorbar(surfplot,shrink=0.75)
#       plt.savefig(ck+"/"+key+"_3D_yz_t_"+tstring+"_old.png",bbox_inches='tight')
          plt.savefig(ck+"/"+key+"_3D_yz_t_"+tstring+".png")
          fig.clf() 
          plt.close()

        if(xzf):
          fig = plt.figure()
          fig.set_size_inches(7, 7)
          ax = fig.gca(projection = '3d')
          if(coord==1):
            ax.set_xlabel('z')
          else:
            ax.set_xlabel(r'$\eta$')
          ax.set_ylabel('x [fm]')
          ax.set_zlabel(ylab)
          X2, Z2 = np.meshgrid(z, x)
          plt.title(key+", t="+tstring+", y="+yc_string)
          plotarr=myarr[:,yy(0),:]
          surfplot=ax.plot_surface(X2, Z2, plotarr, cmap=clm)
          cbar=plt.colorbar(surfplot,shrink=0.75)
#       plt.savefig(ck+"/"+key+"_3D_xz_t_"+tstring+"_old.png",bbox_inches='tight')
          plt.savefig(ck+"/"+key+"_3D_xz_t_"+tstring+".png")
          fig.clf() 
          plt.close()

        if(xpf):
          plt.xlabel('x [fm]')
          plt.ylabel(ylab)
          plt.title(key+", t="+tstring+", y="+yc_string+", z="+zc_string)
          plt.plot(x,myarr[:,yy(0),zz(0)])
          plt.savefig(ck+"/"+key+"_vs_x_t_"+tstring+".png",bbox_inches='tight')
          plt.clf()
          plt.close()

        if(ypf):
          plt.xlabel('y [fm]')
          plt.ylabel(ylab)
          plt.title(key+", t="+tstring+", x="+xc_string+", z="+zc_string)
          plt.plot(y,myarr[xx(0),:,zz(0)])
          plt.savefig(ck+"/"+key+"_vs_y_t_"+tstring+".png",bbox_inches='tight')
          plt.clf()
          plt.close()

        if(zpf):
          if(coord==1):
            plt.xlabel('z')
          else:
            plt.xlabel(r'$\eta$')
          plt.ylabel(ylab)
          plt.title(key+", t="+tstring+", x="+xc_string+", y="+yc_string)
          plt.plot(z,myarr[xx(0),yy(0),:])
          plt.savefig(ck+"/"+key+"_vs_z_t_"+tstring+".png",bbox_inches='tight')
          plt.clf()
          plt.close()
      
        if(xyf):  
          fig = plt.figure(figsize=(7, 7))
          ax = fig.gca()
          plt.gca().set_aspect("equal")
          ax.set_xlabel('x [fm]')
          ax.set_ylabel('y [fm]')
          X2, Y2 = np.meshgrid(x, y)
          plt.title(key+", t="+tstring+", z="+zc_string)
          if(clm==cmap_div):
            parr=myarr[:,:,zz(0)].transpose()
            min_parr=np.amin(parr)
            max_parr=np.amax(parr)
            if((min_parr <= 0) and (max_parr >= 0)):
                if(-min_parr > max_parr):
                    extr_parr=-min_parr
                else:
                    extr_parr=max_parr
#           refval=np.mean(np.abs(myarr[:,:,zz(0)]))
            surfplot=plt.pcolormesh(X2, Y2, myarr[:,:,zz(0)].transpose(), cmap=clm, vmin=-extr_parr, vmax=extr_parr)
          else:
            surfplot=plt.pcolormesh(X2, Y2, myarr[:,:,zz(0)].transpose(), cmap=clm)
          cbar=plt.colorbar(surfplot,shrink=0.67)
          plt.savefig(ck+"/"+key+"_2D_xy_t_"+tstring+".png",bbox_inches='tight')
          fig.clf() 
          plt.close()

        if(yzf):
          fig = plt.figure(figsize=(7, 7))
          ax = fig.gca()
          plt.gca().set_aspect("equal")
          if(coord==1):
            ax.set_ylabel('z')
          else:
            ax.set_ylabel(r'$\eta$')
          ax.set_xlabel('y [fm]')
          Y2, Z2 = np.meshgrid(y, z)
          plt.title(key+", t="+tstring+", x="+xc_string)
          if(clm==cmap_div):
#          refval=np.mean(np.abs(myarr[xx(0):,:]))
            parr=myarr[xx(0),:,:].transpose()
            min_parr=np.amin(parr)
            max_parr=np.amax(parr)
            if((min_parr <= 0) and (max_parr >= 0)):
                if(-min_parr > max_parr):
                    extr_parr=-min_parr
                else:
                    extr_parr=max_parr
            surfplot=plt.pcolormesh(Y2, Z2, myarr[xx(0),:,:].transpose(), cmap=clm, vmin=-extr_parr, vmax=extr_parr)
          else:
            surfplot=plt.pcolormesh(Y2, Z2, myarr[xx(0),:,:].transpose(), cmap=clm)
          cbar=plt.colorbar(surfplot,shrink=0.67)
          plt.savefig(ck+"/"+key+"_2D_yz_t_"+tstring+".png",bbox_inches='tight')
          fig.clf() 
          plt.close()

        if(xzf):
          fig = plt.figure(figsize=(7, 7))
          ax = fig.gca()
          plt.gca().set_aspect("equal")
          if(coord==1):
            ax.set_ylabel('z')
          else:
            ax.set_ylabel(r'$\eta$')
          ax.set_xlabel('x [fm]')
          X2, Z2 = np.meshgrid(x, z)
          plt.title(key+", t="+tstring+", y="+yc_string)
#        plotarr=myarr[:,yy(0),:]
          if(clm==cmap_div):
            parr=myarr[:,yy(0),:].transpose()
            min_parr=np.amin(parr)
            max_parr=np.amax(parr)
            if((min_parr <= 0) and (max_parr >= 0)):
                if(-min_parr > max_parr):
                    extr_parr=-min_parr
                else:
                    extr_parr=max_parr
            surfplot=plt.pcolormesh(X2, Z2, myarr[:,yy(0),:].transpose(), cmap=clm, vmin=-extr_parr, vmax=extr_parr)
          else:
            surfplot=plt.pcolormesh(X2, Z2, myarr[:,yy(0),:].transpose(), cmap=clm)
          cbar=plt.colorbar(surfplot,shrink=0.67)
          plt.savefig(ck+"/"+key+"_2D_xz_t_"+tstring+".png",bbox_inches='tight')
          fig.clf() 
          plt.close()
 
    return 0

def st(c1,c2,v1,v2):
    """It diplays a streamplot of two components of a vector variable, like the velocity or the magnetic field. The streamlines are colored according to the square root of the sum of the component squares (please, note that in general this is not really the module of the variable, because we miss a component and we might miss a factor on the z component when using Milne coordinates.\nSyntax: st(<coord array 1> <coord array 2> <var component 1> <var component 2>)\nInput arguments: the 1D numpy arrays of coordinates of the first and the second component of the variable to be plotted and the variable components given as two 2D numpy arrays.\nReturn value: 0 in case of success, an error will be printed otherwise and/or no plot will be displayed.\nExample: st(x,y,bx[:,:,zz(0)],by[:,:,zz(0)])\n
    """
    import matplotlib.pyplot as plt
    ysn,xsn=np.meshgrid(c1,c2)
    fig = plt.figure(figsize=(6, 6))
    ax1 = fig.add_subplot(1,1,1)
    strm = ax1.streamplot(xsn.T, ysn.T, v1, v2, color=np.sqrt(v1**2+v2**2), linewidth=1, cmap='plasma')
    fig.colorbar(strm.lines,shrink=0.67)
    ax1.set_aspect('equal')
    plt.tight_layout()
    plt.show()
    return 0



def toB(outfile_name, tau, xside, nx_b, yside, ny_b, etaside, neta, chunks=1):
  """It converts the results of a simulation from Minkowski to Bjorken/Milne coordinates with a user chosen grid using the various output files and saves them into an outputfile. The user grid has to be compatible with the Minkowski space time grid used in the simulation, which is reconstructed from the file config_summary.dat.\nSyntax: toB("output file",tau, xside, nx, yside, ny, etaside, neta)\nInput arguments:\n"output file" : the name of the output file"\ntau: the time tau coordinate at which to perform the conversion\nxside: the extension of the grid along positive x in fm (it is assumed that the grid will extend from -xside to +xside\nnx: the number of points of the grid in the x direction\nyside: the same as for xside, but for y\nny: the number of points of the grid in the y direction\netaside: the same as for yside, but for eta\nneta: the number of points of the grid in the eta direction\nReturn value: 0 in case of success, 1 if the output data are already in Milne/Bjorken coordinates.\nExample: dan.toB("Milne_output.dat",1,15,100,15,100,10,100)\n"""

  global x,y,z,rho,vx,vy,vz,en,pr,temp,entr,bulk,pitt,pitx,pity,pitz,pixx,pixy,pixz,piyz,piyy,pizz,g,bx,by,bz,ex,ey,ez,glm,rc
  global t,nx,ny,nz,ncells,xmin,xmax,ymin,ymax,coord,prec,vis,mhd,tfo,xstep,ystep,zstep,var_list

  #we check that we have the essential file config_summary.dat
  try:
    fp=open("config_summary.dat","r")
  except FileNotFoundError:
    print("Sorry, I cannot find the file 'config_summary.dat'.\n")
    sys.exit(1)

  #we store the informations in a string and we close the file  
  confdata=fp.read()
  fp.close()

  #we get informations about the type of simulation
  ss=confdata.find("COORD....=")
  if(ss != -1):
    ee=confdata.find("\n",ss)
    coord=int(confdata[ss+10:ee].strip())
  else:
    print("Error, I cannot determine COORD, i.e. the coordinate type\n")
    sys.exit(1)
  #we check that the info about the coordinate type makes sens
  if(coord==1):
    "Data about simulation in Minkowski coordinates. OK, we can proceed."
  elif(coord==2):
    "Data about simulation alredy in Milne/Bjorken coordinates!!!"
    "No need for conversion."
    return 1
  else:
    "Data about simulation in unknown coordinate system... I quit."
    sys.exit(2)

  #We determine how large is the grid:
  ss=confdata.find("NX.......=")
  if(ss != -1):
    ee=confdata.find("\n",ss)
    nx=int(confdata[ss+10:ee].strip())
  else:
    print("Error, I cannot determine NX\n")
    sys.exit(1)

  ss=confdata.find("NY.......=")
  if(ss != -1):
    ee=confdata.find("\n",ss)
    ny=int(confdata[ss+10:ee].strip())
  else:
    print("Error, I cannot determine NY\n")
    sys.exit(1)

  ss=confdata.find("NZ.......=")
  if(ss != -1):
    ee=confdata.find("\n",ss)
    nz=int(confdata[ss+10:ee].strip())
  else:
    print("Error, I cannot determine NZ\n")
    sys.exit(1)

  ncells=nx*ny*nz

  ss=confdata.find("XMIN.....=")
  if(ss != -1):
    ee=confdata.find("\n",ss)
    xmin=np.float64(confdata[ss+10:ee].strip())
  else:
    print("Error, I cannot determine XMIN\n")
    sys.exit(1)

  ss=confdata.find("XMAX.....=")
  if(ss != -1):
    ee=confdata.find("\n",ss)
    xmax=np.float64(confdata[ss+10:ee].strip())
  else:
    print("Error, I cannot determine XMAX\n")
    sys.exit(1)
  
  ss=confdata.find("YMIN.....=")
  if(ss != -1):
    ee=confdata.find("\n",ss)
    ymin=np.float64(confdata[ss+10:ee].strip())
  else:
    print("Error, I cannot determine YMIN\n")
    sys.exit(1)

  ss=confdata.find("YMAX.....=")
  if(ss != -1):
    ee=confdata.find("\n",ss)
    ymax=np.float64(confdata[ss+10:ee].strip())
  else:
    print("Error, I cannot determine YMAX\n")
    sys.exit(1)

  ss=confdata.find("ZMIN.....=")
  if(ss != -1):
    ee=confdata.find("\n",ss)
    zmin=np.float64(confdata[ss+10:ee].strip())
  else:
    print("Error, I cannot determine ZMIN\n")
    sys.exit(1)

  ss=confdata.find("ZMAX.....=")
  if(ss != -1):
    ee=confdata.find("\n",ss)
    zmax=np.float64(confdata[ss+10:ee].strip())
  else:
    print("Error, I cannot determine ZMAX\n")
    sys.exit(1)

  ss=confdata.find("TSTART...=")
  if(ss != -1):
    ee=confdata.find("\n",ss)
    tstart=np.float64(confdata[ss+10:ee].strip())
  else:
    print("Error, I cannot determine TSTART\n")
    sys.exit(1)

  ss=confdata.find("TSTOP....=")
  if(ss != -1):
    ee=confdata.find("\n",ss)
    tend=np.float64(confdata[ss+10:ee].strip())
  else:
    print("Error, I cannot determine TSTOP\n")
    sys.exit(1)

  ss=confdata.find("DTOUT....=")
  if(ss != -1):
    ee=confdata.find("\n",ss)
    dt_m=np.float64(confdata[ss+10:ee].strip())
  else:
    print("Error, I cannot determine DTOUT\n")
    sys.exit(1)

  ss=confdata.find("OUTP_PREC=")
  if(ss != -1):
    ee=confdata.find("\n",ss)
    prec=confdata[ss+10:ee].strip()
    if(prec == "8"):
      otype=np.float64
      dimfloat=8
      print("Output precision: 8 bytes (double)\n")
    else:
      otype=np.float32
      dimfloat=4
      print("Output precision: 4 bytes (single)\n")
  else:
    print("Error, I cannot determine the output precision (i.e. 8 or 4 bytes)\n")
    sys.exit(1)

  ss=confdata.find("VISCOUS..=")
  if(ss != -1):
    ee=confdata.find("\n",ss)
    vis_flag=confdata[ss+10:ee].strip()
    if(vis_flag == "1"):
      print("Sorry, at the moment this conversion utility works only in the ideal case")
      return 1
    else:
      vis=False
  else:
    print("Error, I cannot determine VISCOUS, i.e. if the simulation is viscous or not\n")
    sys.exit(1)

  ss=confdata.find("MHD......=")
  if(ss != -1):
    ee=confdata.find("\n",ss)
    mhd_flag=confdata[ss+10:ee].strip()
    if(mhd_flag == "1"):
      mhd=True
    else:
      print("Sorry, but at the moment this routine works only in the ideal MHD case\n")
      return 1
  else:
    print("Error, I cannot determine MHD\n")
    sys.exit(1)

  dx_m=(xmax-xmin)/nx
  xgrid_m=np.linspace(xmin+dx_m/2,xmax-dx_m/2,nx)
  dy_m=(ymax-ymin)/ny
  ygrid_m=np.linspace(ymin+dy_m/2,ymax-dy_m/2,ny)
  dz_m=(zmax-zmin)/nz
  zgrid_m=np.linspace(zmin+dz_m/2,zmax-dz_m/2,nz)

  dx_b=2.*xside/nx_b
  xgrid_b=np.linspace(-xside+dx_b/2,xside-dx_b/2,nx_b)
  dy_b=2.*yside/ny_b
  ygrid_b=np.linspace(-yside+dy_b/2,yside-dy_b/2,ny_b)
  deta=2.*etaside/neta
  etagrid=np.linspace(-etaside+deta/2,etaside-deta/2,neta)

  if(max(xgrid_b) > max(xgrid_m)):
    print("Extrapolation error:\n")
    print("Maximum x value of the new grid: "+str(max(xgrid_b))+"\n")
    print("Maximum x value of the old grid: "+str(max(xgrid_m))+"\n")
    return 1
  if(min(xgrid_b) < min(xgrid_m)):
    print("Extrapolation error:\n")
    print("Minimum x value of the new grid: "+str(min(xgrid_b))+"\n")
    print("Minimum x value of the old grid: "+str(min(xgrid_m))+"\n")
    return 1

  if(max(ygrid_b) > max(ygrid_m)):
    print("Extrapolation error:\n")
    print("Maximum y value of the new grid: "+str(max(ygrid_b))+"\n")
    print("Maximum y value of the old grid: "+str(max(ygrid_m))+"\n")
    return 1
  if(min(ygrid_b) < min(ygrid_m)):
    print("Extrapolation error:\n")
    print("Minimum y value of the new grid: "+str(min(ygrid_b))+"\n")
    print("Minimum y value of the old grid: "+str(min(ygrid_m))+"\n")
    return 1

  if(tau*math.sinh(max(etagrid)) > max(zgrid_m)):
    print("Extrapolation error:\n")
    print("Maximum corresponding z value of the new grid: "+str(tau*math.sinh(max(etagrid)))+"\n")
    print("Maximum z value of the old grid: "+str(max(zgrid_m))+"\n")
    return 1
  if(tau*math.sinh(min(etagrid)) < min(zgrid_m)):
    print("Extrapolation error:\n")
    print("Maximum corresponding z value of the new grid: "+str(tau*math.sinh(min(etagrid)))+"\n")
    print("Minimum z value of the old grid: "+str(min(zgrid_m))+"\n")
    return 1

  if(tau*math.cosh(max(etagrid)) > tend):
    #we do not need to check the minimum as etagrid has been built symmetric
    print("Extrapolation error:\n")
    print("Maximum corresponding t value of the new grid: "+str(tau*math.cosh(max(etagrid)))+"\n")
    print("Maximum t value of the old grid: "+str(tend)+"\n")
    return 1

  ttop=tau*math.cosh(max(etagrid))

  nt=int(math.ceil((ttop-tstart)/dt_m))+1

  tgrid=np.linspace(tstart,tstart+(nt-1)*dt_m,num=nt)

  ut_in=np.zeros((nt,nx,ny,nz),dtype=np.float64)
  ux_in=np.zeros((nt,nx,ny,nz),dtype=np.float64)
  uy_in=np.zeros((nt,nx,ny,nz),dtype=np.float64)
  uz_in=np.zeros((nt,nx,ny,nz),dtype=np.float64)
  pr_in=np.zeros((nt,nx,ny,nz),dtype=np.float64)
  bx_in=np.zeros((nt,nx,ny,nz),dtype=np.float64)
  by_in=np.zeros((nt,nx,ny,nz),dtype=np.float64)
  bz_in=np.zeros((nt,nx,ny,nz),dtype=np.float64)
  ex_in=np.zeros((nt,nx,ny,nz),dtype=np.float64)
  ey_in=np.zeros((nt,nx,ny,nz),dtype=np.float64)
  ez_in=np.zeros((nt,nx,ny,nz),dtype=np.float64)
 
  for i in range(0,nt):
      infile="out"+'{:04d}'.format(i+1)+".dat"
      load(infile,"vx","vy","vz","pr","bx","by","bz","ex","ey","ez")
      glf=1./np.sqrt(1-vx**2-vy**2-vz**2)
      ut_in[i,:,:,:]=glf
      ux_in[i,:,:,:]=vx*glf
      uy_in[i,:,:,:]=vy*glf
      uz_in[i,:,:,:]=vz*glf
      pr_in[i,:,:,:]=pr
      bx_in[i,:,:,:]=bx
      by_in[i,:,:,:]=by
      bz_in[i,:,:,:]=bz
      ex_in[i,:,:,:]=ex
      ey_in[i,:,:,:]=ey
      ez_in[i,:,:,:]=ez


  ut_interpolator = intp.RegularGridInterpolator((tgrid, xgrid_m, ygrid_m, zgrid_m), ut_in, fill_value=0.)
  ux_interpolator = intp.RegularGridInterpolator((tgrid, xgrid_m, ygrid_m, zgrid_m), ux_in, fill_value=0.)
  uy_interpolator = intp.RegularGridInterpolator((tgrid, xgrid_m, ygrid_m, zgrid_m), uy_in, fill_value=0.)
  uz_interpolator = intp.RegularGridInterpolator((tgrid, xgrid_m, ygrid_m, zgrid_m), uz_in, fill_value=0.)
  pr_interpolator = intp.RegularGridInterpolator((tgrid, xgrid_m, ygrid_m, zgrid_m), pr_in, fill_value=0.)
  bx_interpolator = intp.RegularGridInterpolator((tgrid, xgrid_m, ygrid_m, zgrid_m), bx_in, fill_value=0.)
  by_interpolator = intp.RegularGridInterpolator((tgrid, xgrid_m, ygrid_m, zgrid_m), by_in, fill_value=0.)
  bz_interpolator = intp.RegularGridInterpolator((tgrid, xgrid_m, ygrid_m, zgrid_m), bz_in, fill_value=0.)
  ex_interpolator = intp.RegularGridInterpolator((tgrid, xgrid_m, ygrid_m, zgrid_m), ex_in, fill_value=0.)
  ey_interpolator = intp.RegularGridInterpolator((tgrid, xgrid_m, ygrid_m, zgrid_m), ey_in, fill_value=0.)
  ez_interpolator = intp.RegularGridInterpolator((tgrid, xgrid_m, ygrid_m, zgrid_m), ez_in, fill_value=0.)

  #we allocate the output arrays
  vx_b=np.zeros((nx_b,ny_b,neta),dtype=np.float64)
  vy_b=np.zeros((nx_b,ny_b,neta),dtype=np.float64)
  vz_b=np.zeros((nx_b,ny_b,neta),dtype=np.float64)
  pr_b=np.zeros((nx_b,ny_b,neta),dtype=np.float64)
  bx_b=np.zeros((nx_b,ny_b,neta),dtype=np.float64)
  by_b=np.zeros((nx_b,ny_b,neta),dtype=np.float64)
  bz_b=np.zeros((nx_b,ny_b,neta),dtype=np.float64)
  ex_b=np.zeros((nx_b,ny_b,neta),dtype=np.float64)
  ey_b=np.zeros((nx_b,ny_b,neta),dtype=np.float64)
  ez_b=np.zeros((nx_b,ny_b,neta),dtype=np.float64)

  outfile=open(outfile_name,"w")
  sp="    "
  for i in range(0,nx_b):
      x_b=xgrid_b[i]
      for j in range(0,ny_b):
          y_b=ygrid_b[j]
          for k in range(0,neta):
              eta=etagrid[k]
              ch=math.cosh(eta)
              sh=math.sinh(eta)
              tt=tau*ch
              zz=tau*sh
              ut_m=ut_interpolator((tt,x_b,y_b,zz))
              ux_m=ux_interpolator((tt,x_b,y_b,zz))
              uy_m=uy_interpolator((tt,x_b,y_b,zz))
              uz_m=uz_interpolator((tt,x_b,y_b,zz))
              pr_m=pr_interpolator((tt,x_b,y_b,zz))
              bx_m=bx_interpolator((tt,x_b,y_b,zz))
              by_m=by_interpolator((tt,x_b,y_b,zz))
              bz_m=bz_interpolator((tt,x_b,y_b,zz))
              ex_m=ex_interpolator((tt,x_b,y_b,zz))
              ey_m=ey_interpolator((tt,x_b,y_b,zz))
              ez_m=ez_interpolator((tt,x_b,y_b,zz))
              ut_b=ut_m*ch-uz_m*sh
              ux_b=ux_m #just for clarity
              uy_b=uy_m #just for clarity
              uz_b=(-ut_m*sh+uz_m*ch)/tau
              vx_b[i,j,k]=ux_b/ut_b
              vy_b[i,j,k]=uy_b/ut_b
              vz_b[i,j,k]=uz_b/ut_b
              pr_b[i,j,k]=pr_m
              bx_b[i,j,k]=bx_m*ch+ey_m*sh
              by_b[i,j,k]=by_m*ch-ex_m*sh
              bz_b[i,j,k]=bz_m/tau
              ex_b[i,j,k]=ex_m*ch-by_m*sh
              ey_b[i,j,k]=ey_m*ch+bx_m*sh
              ez_b[i,j,k]=ex_m/tau
              outfile.write(str(x_b)+sp+str(y_b)+sp+str(eta)+sp+str(vx_b[i,j,k])+sp+str(vy_b[i,j,k])+sp+str(vz_b[i,j,k])+sp+str(pr_b[i,j,k])+sp+str(bx_b[i,j,k])+sp+str(by_b[i,j,k])+sp+str(bz_b[i,j,k])+"\n")
  
  outfile.close()
  print(outfile_name+" written.")                      


def toM(outfile_name, t_m, xside, nx_m, yside, ny_m, zside, nz_m, chunks=1):
  """It converts the results of a simulation from Bjorken/Milne to Minkowski coordinates with a user chosen grid using the various output files and saves them into an outputfile. The user grid has to be compatible with the Bjorken/Milne space time grid used in the simulation, which is reconstructed from the file config_summary.dat.\nSyntax: toM("output file",t, xside, nx, yside, ny, zside, nz)\nInput arguments:\n"output file" : the name of the output file"\nt: the time coordinate at which to perform the conversion\nxside: the extension of the grid along positive x in fm (it is assumed that the grid will extend from -xside to +xside\nnx: the number of points of the grid in the x direction\nyside: the same as for xside, but for y\nny: the number of points of the grid in the y direction\nzside: the same as for yside, but for z\nnz: the number of points of the grid in the z direction\nReturn value: 0 in case of success, 1 if the output data are already in Minkowski coordinates.\nExample: dan.toM("Minkowski_output.dat",1,15,100,15,100,10,100)\n"""

  global x,y,z,rho,vx,vy,vz,en,pr,temp,entr,bulk,pitt,pitx,pity,pitz,pixx,pixy,pixz,piyz,piyy,pizz,g,bx,by,bz,ex,ey,ez,glm,rc
  global t,nx,ny,nz,ncells,xmin,xmax,ymin,ymax,coord,prec,vis,mhd,tfo,xstep,ystep,zstep,var_list

  #we check that we have the essential file config_summary.dat
  try:
    fp=open("config_summary.dat","r")
  except FileNotFoundError:
    print("Sorry, I cannot find the file 'config_summary.dat'.\n")
    sys.exit(1)

  #we store the informations in a string and we close the file  
  confdata=fp.read()
  fp.close()

  #we get informations about the type of simulation
  ss=confdata.find("COORD....=")
  if(ss != -1):
    ee=confdata.find("\n",ss)
    coord=int(confdata[ss+10:ee].strip())
  else:
    print("Error, I cannot determine COORD, i.e. the coordinate type\n")
    sys.exit(1)
  #we check that the info about the coordinate type makes sens
  if(coord==2):
    "Data about simulation in Milne/Bjorken coordinates. OK, we can proceed."
  elif(coord==2):
    "Data about simulation alredy in Minkowski coordinates!!!"
    "No need for conversion."
    return 1
  else:
    "Data about simulation in unknown coordinate system... I quit."
    sys.exit(2)

  #We determine how large is the grid:
  ss=confdata.find("NX.......=")
  if(ss != -1):
    ee=confdata.find("\n",ss)
    nx_b=int(confdata[ss+10:ee].strip())
  else:
    print("Error, I cannot determine NX\n")
    sys.exit(1)

  ss=confdata.find("NY.......=")
  if(ss != -1):
    ee=confdata.find("\n",ss)
    ny_b=int(confdata[ss+10:ee].strip())
  else:
    print("Error, I cannot determine NY\n")
    sys.exit(1)

  ss=confdata.find("NZ.......=")
  if(ss != -1):
    ee=confdata.find("\n",ss)
    nz_b=int(confdata[ss+10:ee].strip())
  else:
    print("Error, I cannot determine NZ\n")
    sys.exit(1)

  ss=confdata.find("XMIN.....=")
  if(ss != -1):
    ee=confdata.find("\n",ss)
    xmin_b=np.float64(confdata[ss+10:ee].strip())
  else:
    print("Error, I cannot determine XMIN\n")
    sys.exit(1)

  ss=confdata.find("XMAX.....=")
  if(ss != -1):
    ee=confdata.find("\n",ss)
    xmax_b=np.float64(confdata[ss+10:ee].strip())
  else:
    print("Error, I cannot determine XMAX\n")
    sys.exit(1)
  
  ss=confdata.find("YMIN.....=")
  if(ss != -1):
    ee=confdata.find("\n",ss)
    ymin_b=np.float64(confdata[ss+10:ee].strip())
  else:
    print("Error, I cannot determine YMIN\n")
    sys.exit(1)

  ss=confdata.find("YMAX.....=")
  if(ss != -1):
    ee=confdata.find("\n",ss)
    ymax_b=np.float64(confdata[ss+10:ee].strip())
  else:
    print("Error, I cannot determine YMAX\n")
    sys.exit(1)

  ss=confdata.find("ZMIN.....=")
  if(ss != -1):
    ee=confdata.find("\n",ss)
    zmin_b=np.float64(confdata[ss+10:ee].strip())
  else:
    print("Error, I cannot determine ZMIN\n")
    sys.exit(1)

  ss=confdata.find("ZMAX.....=")
  if(ss != -1):
    ee=confdata.find("\n",ss)
    zmax_b=np.float64(confdata[ss+10:ee].strip())
  else:
    print("Error, I cannot determine ZMAX\n")
    sys.exit(1)

  ss=confdata.find("TSTART...=")
  if(ss != -1):
    ee=confdata.find("\n",ss)
    tstart=np.float64(confdata[ss+10:ee].strip())
  else:
    print("Error, I cannot determine TSTART\n")
    sys.exit(1)

  ss=confdata.find("TSTOP....=")
  if(ss != -1):
    ee=confdata.find("\n",ss)
    tend=np.float64(confdata[ss+10:ee].strip())
  else:
    print("Error, I cannot determine TSTOP\n")
    sys.exit(1)

  ss=confdata.find("DTOUT....=")
  if(ss != -1):
    ee=confdata.find("\n",ss)
    dt_b=np.float64(confdata[ss+10:ee].strip())
  else:
    print("Error, I cannot determine DTOUT\n")
    sys.exit(1)

  ss=confdata.find("OUTP_PREC=")
  if(ss != -1):
    ee=confdata.find("\n",ss)
    prec=confdata[ss+10:ee].strip()
    if(prec == "8"):
      otype=np.float64
      dimfloat=8
      print("Output precision: 8 bytes (double)\n")
    else:
      otype=np.float32
      dimfloat=4
      print("Output precision: 4 bytes (single)\n")
  else:
    print("Error, I cannot determine the output precision (i.e. 8 or 4 bytes)\n")
    sys.exit(1)

  ss=confdata.find("VISCOUS..=")
  if(ss != -1):
    ee=confdata.find("\n",ss)
    vis_flag=confdata[ss+10:ee].strip()
    if(vis_flag == "1"):
      print("Sorry, at the moment this conversion utility works only in the ideal case")
      return 1
    else:
      vis=False
  else:
    print("Error, I cannot determine VISCOUS, i.e. if the simulation is viscous or not\n")
    sys.exit(1)

  ss=confdata.find("MHD......=")
  if(ss != -1):
    ee=confdata.find("\n",ss)
    mhd_flag=confdata[ss+10:ee].strip()
    if(mhd_flag == "1"):
      mhd=True
    else:
      print("Sorry, but at the moment this routine works only in the ideal MHD case\n")
      return 1
  else:
    print("Error, I cannot determine MHD\n")
    sys.exit(1)

  dx_b=(xmax_b-xmin_b)/nx_b
  xgrid_b=np.linspace(xmin_b+dx_b/2,xmax_b-dx_b/2,nx_b)
  dy_b=(ymax_b-ymin_b)/ny_b
  ygrid_b=np.linspace(ymin_b+dy_b/2,ymax_b-dy_b/2,ny_b)
  dz_b=(zmax_b-zmin_b)/nz_b
  zgrid_b=np.linspace(zmin_b+dz_b/2,zmax_b-dz_b/2,nz_b)

  dx_m=2.*xside/nx_m
  xgrid_m=np.linspace(-xside+dx_m/2,xside-dx_m/2,nx_m)
  dy_m=2.*yside/ny_m
  ygrid_m=np.linspace(-yside+dy_m/2,yside-dy_m/2,ny_m)
  dz_m=2.*zside/nz_m
  zgrid_m=np.linspace(-zside+dz_m/2,zside-dz_m/2,nz_m)

  if(max(xgrid_m) > max(xgrid_b)):
    print("Extrapolation error:\n")
    print("Maximum x value of the new grid: "+str(max(xgrid_m))+"\n")
    print("Maximum x value of the old grid: "+str(max(xgrid_b))+"\n")
    return 1
  if(min(xgrid_m) < min(xgrid_b)):
    print("Extrapolation error:\n")
    print("Minimum x value of the new grid: "+str(min(xgrid_m))+"\n")
    print("Minimum x value of the old grid: "+str(min(xgrid_b))+"\n")
    return 1

  if(max(ygrid_m) > max(ygrid_b)):
    print("Extrapolation error:\n")
    print("Maximum y value of the new grid: "+str(max(ygrid_m))+"\n")
    print("Maximum y value of the old grid: "+str(max(ygrid_b))+"\n")
    return 1
  if(min(ygrid_m) < min(ygrid_b)):
    print("Extrapolation error:\n")
    print("Minimum y value of the new grid: "+str(min(ygrid_m))+"\n")
    print("Minimum y value of the old grid: "+str(min(ygrid_b))+"\n")
    return 1

  if(max(zgrid_m) > max(zgrid_b)):
    print("Extrapolation error:\n")
    print("Maximum corresponding z value of the new grid: "+str(max(zgrid_m))+"\n")
    print("Maximum z value of the old grid: "+str(max(zgrid_b))+"\n")
    return 1
  if(min(zgrid_m) < min(zgrid_b)):
    print("Extrapolation error:\n")
    print("Maximum corresponding z value of the new grid: "+str(min(zgrid_m))+"\n")
    print("Minimum z value of the old grid: "+str(min(zgrid_b))+"\n")
    return 1

  if(t_m > tend):
    #we do not need to check the minimum as etagrid has been built symmetric
    print("Extrapolation error:\n")
    print("Desired time in Minkowski coordinates "+str(t_m)+" is larger than the final time of the simulation in Bjorken coord.: "+str(tend))
    return 1

  nt_b=int(math.ceil((tend-tstart)/dt_b))+1

  print("nt_b is:",nt_b)

  tgrid=np.linspace(tstart,tstart+(nt_b-1)*dt_b,num=nt_b)

  print("tstart:"+str(tstart)+", tend:"+str(tend)+", dt_b:"+str(dt_b))
  print("tgrid:"+str(tgrid))
  print("xgrid_b:"+str(xgrid_b))
  print("ygrid_b:"+str(ygrid_b))
  print("zgrid_b:"+str(zgrid_b))

  ut_in=np.zeros((nt_b,nx_b,ny_b,nz_b),dtype=np.float64)
  ux_in=np.zeros((nt_b,nx_b,ny_b,nz_b),dtype=np.float64)
  uy_in=np.zeros((nt_b,nx_b,ny_b,nz_b),dtype=np.float64)
  uz_in=np.zeros((nt_b,nx_b,ny_b,nz_b),dtype=np.float64)
  pr_in=np.zeros((nt_b,nx_b,ny_b,nz_b),dtype=np.float64)
  bx_in=np.zeros((nt_b,nx_b,ny_b,nz_b),dtype=np.float64)
  by_in=np.zeros((nt_b,nx_b,ny_b,nz_b),dtype=np.float64)
  bz_in=np.zeros((nt_b,nx_b,ny_b,nz_b),dtype=np.float64)
  ex_in=np.zeros((nt_b,nx_b,ny_b,nz_b),dtype=np.float64)
  ey_in=np.zeros((nt_b,nx_b,ny_b,nz_b),dtype=np.float64)
  ez_in=np.zeros((nt_b,nx_b,ny_b,nz_b),dtype=np.float64)
 
  for i in range(0,nt_b):
      infile="out"+'{:04d}'.format(i+1)+".dat"
      load(infile,"vx","vy","vz","pr","bx","by","bz","ex","ey","ez")
      glf=1./np.sqrt(1-vx**2-vy**2-(vz*t)**2)
      ut_in[i,:,:,:]=glf
      ux_in[i,:,:,:]=vx*glf
      uy_in[i,:,:,:]=vy*glf
      uz_in[i,:,:,:]=vz*glf
      pr_in[i,:,:,:]=pr
      bx_in[i,:,:,:]=bx
      by_in[i,:,:,:]=by
      bz_in[i,:,:,:]=bz
      ex_in[i,:,:,:]=ex
      ey_in[i,:,:,:]=ey
      ez_in[i,:,:,:]=ez


  ut_interpolator = intp.RegularGridInterpolator((tgrid, xgrid_b, ygrid_b, zgrid_b), ut_in, fill_value=0.)
  ux_interpolator = intp.RegularGridInterpolator((tgrid, xgrid_b, ygrid_b, zgrid_b), ux_in, fill_value=0.)
  uy_interpolator = intp.RegularGridInterpolator((tgrid, xgrid_b, ygrid_b, zgrid_b), uy_in, fill_value=0.)
  uz_interpolator = intp.RegularGridInterpolator((tgrid, xgrid_b, ygrid_b, zgrid_b), uz_in, fill_value=0.)
  pr_interpolator = intp.RegularGridInterpolator((tgrid, xgrid_b, ygrid_b, zgrid_b), pr_in, fill_value=0.)
  bx_interpolator = intp.RegularGridInterpolator((tgrid, xgrid_b, ygrid_b, zgrid_b), bx_in, fill_value=0.)
  by_interpolator = intp.RegularGridInterpolator((tgrid, xgrid_b, ygrid_b, zgrid_b), by_in, fill_value=0.)
  bz_interpolator = intp.RegularGridInterpolator((tgrid, xgrid_b, ygrid_b, zgrid_b), bz_in, fill_value=0.)
  ex_interpolator = intp.RegularGridInterpolator((tgrid, xgrid_b, ygrid_b, zgrid_b), ex_in, fill_value=0.)
  ey_interpolator = intp.RegularGridInterpolator((tgrid, xgrid_b, ygrid_b, zgrid_b), ey_in, fill_value=0.)
  ez_interpolator = intp.RegularGridInterpolator((tgrid, xgrid_b, ygrid_b, zgrid_b), ez_in, fill_value=0.)

  #we allocate the output arrays
  vx_m=np.zeros((nx_m,ny_m,nz_m),dtype=np.float64)
  vy_m=np.zeros((nx_m,ny_m,nz_m),dtype=np.float64)
  vz_m=np.zeros((nx_m,ny_m,nz_m),dtype=np.float64)
  pr_m=np.zeros((nx_m,ny_m,nz_m),dtype=np.float64)
  bx_m=np.zeros((nx_m,ny_m,nz_m),dtype=np.float64)
  by_m=np.zeros((nx_m,ny_m,nz_m),dtype=np.float64)
  bz_m=np.zeros((nx_m,ny_m,nz_m),dtype=np.float64)
  ex_m=np.zeros((nx_m,ny_m,nz_m),dtype=np.float64)
  ey_m=np.zeros((nx_m,ny_m,nz_m),dtype=np.float64)
  ez_m=np.zeros((nx_m,ny_m,nz_m),dtype=np.float64)

  outfile=open(outfile_name,"w")
  sp="    "
  for i in range(0,nx_m):
      x=xgrid_m[i]
      for j in range(0,ny_m):
          y=ygrid_m[j]
          for k in range(0,nz_m):
              z=zgrid_m[k]
              radical=t_m**2-z**2
              if(radical <= 0):
                continue
                #the other condition for eta=1/2 ln((t+z)/(t-z)) is already included here
              tau=math.sqrt(radical)
              eta=0.5*math.log((t_m+z)/(t_m-z))
              ch=math.cosh(eta)
              sh=math.sinh(eta)
#              print("Interpolating: "+str(tau)+sp+str(eta))
              ut_b=ut_interpolator((tau,x,y,eta))
              ux_b=ux_interpolator((tau,x,y,eta))
              uy_b=uy_interpolator((tau,x,y,eta))
              uz_b=uz_interpolator((tau,x,y,eta))
              pr_b=pr_interpolator((tau,x,y,eta))
              bx_b=bx_interpolator((tau,x,y,eta))
              by_b=by_interpolator((tau,x,y,eta))
              bz_b=bz_interpolator((tau,x,y,eta))
              ex_b=ex_interpolator((tau,x,y,eta))
              ey_b=ey_interpolator((tau,x,y,eta))
              ez_b=ez_interpolator((tau,x,y,eta))
              ut_m=ut_b*ch+uz_b*sh*tau
              ux_m=ux_b #just for clarity
              uy_m=uy_b #just for clarity
              uz_m=ut_b*sh+uz_b*ch*tau
              vx_m[i,j,k]=ux_m/ut_m
              vy_m[i,j,k]=uy_m/ut_m
              vz_m[i,j,k]=uz_m/ut_m
              pr_m[i,j,k]=pr_b
              bx_m[i,j,k]=bx_b*ch-ey_b*sh
              by_m[i,j,k]=by_b*ch+ex_b*sh
              bz_m[i,j,k]=bz_b*tau
              ex_m[i,j,k]=ex_b*ch+by_b*sh
              ey_m[i,j,k]=ey_b*ch-bx_b*sh
              ez_m[i,j,k]=ex_b*tau
              outfile.write(str(x)+sp+str(y)+sp+str(z)+sp+str(vx_m[i,j,k])+sp+str(vy_m[i,j,k])+sp+str(vz_m[i,j,k])+sp+str(pr_m[i,j,k])+sp+str(bx_m[i,j,k])+sp+str(by_m[i,j,k])+sp+str(bz_m[i,j,k])+"\n")
  
  outfile.close()
  print(outfile_name+" written.")                      

  outfile=open(outfile_name+"_4gp","w")
  sp="    "
  for k in range(0,nz_m):
      z=zgrid_m[k]
      for i in range(0,nx_m):
          x=xgrid_m[i]
          for j in range(0,ny_m):
              y=ygrid_m[j]
              outfile.write(str(x)+sp+str(y)+sp+str(z)+sp+str(vx_m[i,j,k])+sp+str(vy_m[i,j,k])+sp+str(vz_m[i,j,k])+sp+str(pr_m[i,j,k])+sp+str(bx_m[i,j,k])+sp+str(by_m[i,j,k])+sp+str(bz_m[i,j,k])+"\n")
          outfile.write("\n")
      outfile.write("\n")

  outfile.close()
  print(outfile_name+"_4gp"+" written.")                      

def to_ascii_gp(infile_name,outfile_name):
    """It transforms a binary output file into an ascii file, with spacings suitable to be plotted with gnuplot.\nSyntax: toascii(<name of the ECHO-QGP binary output file to be converted>,<name of the output file>\nExample: toascii("out0001.dat","out1ascii.dat")\n"""
   
    load(infile_name) 
    sp="    "    

    fout=open(outfile_name,"w+")
    fout.write("#  1 x  2 y  3 z")
    for i in range(len(var_list)):
        fout.write("  "+'{:2d}'.format(i+4)+"  "+var_list[i])
    fout.write("\n")
    for iz in range(nz):  
        zp='{:9.6f}'.format(z[iz])
        for ix in range(nx):  
            xp='{:9.6f}'.format(x[ix])+sp
            for iy in range(ny):
                yp='{:9.6f}'.format(y[iy])+sp
                fout.write(xp+yp+zp)
                for il in var_list:
                    myarr=find_arr(il)
                    fout.write(sp+'{:16.14e}'.format(myarr[ix,iy,iz]))
                fout.write("\n")
            fout.write("\n")
        fout.write("\n")

    fout.close()
    return 0 

def check_ideal_gub(output_number, outfile_name):
    """It prints into an ascii output file the values of the temperature, the energy density and the velocity along the diagonal x=y, from 0 to the border of the grid.\nWe assume that the output files follow the default conventions and are named out0001.dat, out0002.dat... and der0001.dat...\nSyntax: check_ideal_gub(<number of the output file to be converted>,<name of the output file>\nExample: check_ideal_gub(1,"gubser1.dat")\n"""

    infile_name="out"+'{:04d}'.format(output_number)+".dat"
    derfile_name="der"+'{:04d}'.format(output_number)+".dat"
    
    print(infile_name)
    print(derfile_name)
    
    
    if(os.path.exists(infile_name)):
        load(infile_name)
    else:
        print(infile_name+" not found")
        return False
	
    if(os.path.exists(derfile_name)):
        load_der(derfile_name)
    else:
        print(derfile_name+" not found")
        return False
    
    sp="    "
    if(nz>1):
        print("Error, I expected to analize a 2D ECHO-QGP output file, but this is 3D. Something went wrong, please, check...")
        print("File "+otputfile_name+" not printed.\n")
        return 0
    if(nx != ny):
        print("Error, I expected an equal number of cells for both nx and ny. Please, check...")
        print("File "+otputfile_name+" not printed.\n")
        return 0
    
    vr=np.sqrt(vx**2+vy**2)
    ur=1./np.sqrt(1-vr**2)
        
    fout=open(outfile_name,"w+")
    fout.write("#  1 x=y   2 en density [MeV/fm^3]   3 temperature [MeV]   4 radial velocity [c units]    5 expansion rate\n")
    for ix in range(xx(0),nx):  
        theta=derivs[dtt,ix,ix,0]+derivs[dxx,ix,ix,0]+derivs[dyy,ix,ix,0]+derivs[dzz,ix,ix,0]+ur[ix,ix,0]/t
        fout.write('{:9.6f}'.format(x[ix]*math.sqrt(2.))+sp+'{:16.14e}'.format(1000*en[ix,ix,0])+sp+'{:16.14e}'.format(1000*temp[ix,ix,0])+sp+'{:16.14e}'.format(vr[ix,ix,0])+sp+'{:16.14e}'.format(theta)+"\n")
    fout.close()
    return True

#now a function to get help
def help():
    print("List of available functions.\nThe arguments in < > brackets are mandatory, the arguments in { } brackets are optional.\nThe list of available variable arrays is in the var_list array, x, y and z are the coordinate 1D arrays and time is stored in the variable t.\n")
    print("xx :  "+xx.__doc__+"\n")
    print("yy :  "+yy.__doc__+"\n")
    print("zz :  "+zz.__doc__+"\n")
    print("load : "+load.__doc__+"\n")
    print("stats :  "+stats.__doc__+"\n") 
    print("print_to_file :  "+print_to_file.__doc__+"\n")
    print("print_for_HF :  "+print_for_HF.__doc__+"\n")
    print("comp_en_hd : "+comp_en_hd.__doc__+"\n")
    print("p3d : "+p3d.__doc__+"\n")
    print("st : "+st.__doc__+"\n")
    print("toB : "+toB.__doc__+"\n")
    print("toM : "+toM.__doc__+"\n")
    print("to_ascii_gp :"+to_ascii_gp.__doc__+"\n")
    print("check_ideal_gub : "+check_ideal_gub.__doc__+"\n")
