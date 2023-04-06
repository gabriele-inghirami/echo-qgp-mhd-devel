# This script makes 2D plots of the initial magnetic and electric
# fields used to initialize ECHO-QGP.
# The plots have common overall maximum and minimum values,
# so they are suitable to make an animation
set term png enh font "Helvetica, 14" size 800,600

set xlabel "x (fm)"
set ylabel "y (fm)"
set cblabel "units of m^2_{{/Symbol p}_0}"

inputfile="EM_z0_4gp"

index_min=0
index_max=0

st=0.0 #the dz step

xside=20
yside=20
zside=1. #the left edge of the z axis 

iH1x=4
iH2x=5
iHx=6
iH1Cx=7
iH2Cx=8
iHCx=9
iHxtot=10
iH1y=11
iH2y=12
iHy=13
iH1Cy=14
iH2Cy=15
iHCy=16
iHytot=17
iH1Cz=18
iH2Cz=19
iHCz=20
iE1x=21
iE2x=22
iEx=23
iE1Cx=24
iE2Cx=25
iECx=26
iExtot=27
iE1y=28
iE2y=29
iEy=30
iE1Cy=31
iE2Cy=32
iECy=33
iEytot=34
iE1z=35
iE2z=36
iEz=37


do for[i=index_min:index_max] {


#pos=-zside+(i+0.5)*st
pos=0

posstring=sprintf( "-z%+-5.2f", pos )
postitlestring=sprintf( "-z=%+-5.2f", pos )

set size square 
set pm3d
unset surface  # don't need surfaces
set view map
#set contour
set key outside
set cntrparam cubicspline  # smooth out the lines
set pm3d interpolate 0,0# interpolate the color
set xrange [-xside:xside]
set yrange [-yside:yside]

set palette defined (-0.05 "#2D4EEE", -0.04 "#4668F0", -0.03  "#809AF4", -0.02 "#C6D4FA", -0.01 "#F2F5F9", 0.01 "#FEFDC6", 0.02 "#FBF463", 0.03 "#EAB540", 0.04 "#D95E31", 0.05 "#D0022A")

set cbrange [-0.05:0.05]

set out "Ext_comp_".posstring.".png"; set title "E_x".postitlestring; splot inputfile index i using 1:2:iExtot notitle with lines lt 1 lw 0; 

set out "Eyt_comp_".posstring.".png"; set title "E_y".postitlestring; splot inputfile index i using 1:2:iEytot notitle with lines lt 1 lw 0; 

set out "Ezt_comp_".posstring.".png"; set title "10 E_z".postitlestring; splot inputfile index i using 1:2:($37)*10 notitle with lines lt 1 lw 0; 


set cbrange [-0.0125:0.0125]

set palette defined (-0.0125 "#2D4EEE", -0.01 "#4668F0", -0.0075 "#809AF4", -0.005 "#C6D4FA", -0.0025 "#F2F5F9", 0.0025 "#FEFDC6", 0.005 "#FBF463", 0.0075 "#EAB540", 0.01 "#D95E31", 0.0125 "#D0022A")

set out "Bxt_comp_".posstring.".png"; set title "B_x".postitlestring; splot inputfile index i using 1:2:iHxtot notitle with lines lt 1 lw 0; 

set out "Byt_comp_".posstring.".png"; set title "B_y/2".postitlestring; splot inputfile index i using 1:2:($17)/2. notitle with lines lt 1 lw 0; 

set out "Bzt_comp_".posstring.".png"; set title "10 B_z".postitlestring; splot inputfile index i using 1:2:($20)*10 notitle with lines lt 1 lw 0; 

}
