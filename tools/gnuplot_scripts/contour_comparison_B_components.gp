# This script makes 2D plots of the x and y components of the magnetic field
# for a certain interval of time

# Set variables
echo_qgp_dir="../../"
dir1=echo_qgp_dir."outr0001/postproc/readx/"
title1="B_x"
title2="B_y"


# Set terminal
set terminal pngcairo size 1400, 700 enhanced font 'Helvetica,22'
set encoding utf8

# Set plot properties
set size square
set pm3d
unset surface 
set view map
set key outside
set pm3d interpolate 0,0# interpolate the color


# Set variable to plot and range of frames
var1="BX2d"
var2="BY2d"
iniframe=1
endframe=201

# Set the axes
set xlabel 'x (fm)'
#set ylabel "{/Symbol h}"
set ylabel 'y (fm)'
set xrange [-15:15]
set yrange [-15:15]


# Set the color palette
# example of method to get stat info directly from datafile
#i0 = sprintf(dir1.var.'%04.0f.dat',iniframe)
#stats i0 using 3
#stmin=STATS_min
#stmax=STATS_max
#step=stmax-stmin

stmin1=system("head -1 " . var1."_stats" . " | awk '{print $1}'")
stmax1=system("head -1 " . var1."_stats" . " | awk '{print $2}'")
stmin2=system("head -1 " . var2."_stats" . " | awk '{print $1}'")
stmax2=system("head -1 " . var2."_stats" . " | awk '{print $2}'")

step1=stmax1-stmin1
step2=stmax2-stmin2

#set palette model RGB defined (stmin 'black', stmin+step/6 'blue', stmin+step/3 'cyan' , stmin+step/2 'green' , stmin+step*2/3. 'yellow', stmin+step*5/6. 'red', stmax 'purple')

do for [i=iniframe:endframe] {
   keytitle="z=0 - t=".sprintf("%5.2f",(i-1)*0.02+0.5)." fm/c"
   infile1 = sprintf(dir1.var1.'%04.0f.dat',i)
   infile2 = sprintf(dir1.var2.'%04.0f.dat',i)
   outfile = sprintf("B_components".'%04.0f.png',i)
   set output outfile
   set multiplot layout 1,2 title keytitle
   set title title1
   set cblabel "B_x [GeV^{1/2}/fm^{3/2}]"
   set cbrange[stmin1:stmax1]
   set palette model RGB defined (stmin1 'black', stmin1+step1/2 'red' , stmax1 'yellow')
   splot infile1 u 1:2:3 notitle
   set title title2
   set cblabel "B_y [GeV^{1/2}/fm^{3/2}]"
   set cbrange[stmin2:stmax2]
   set palette model RGB defined (stmin2 'black', stmin2+step2/2 'red' , stmax2 'yellow')
   splot infile2 u 1:2:3 notitle
   unset multiplot
   print "Frame ".i." done"
}                

