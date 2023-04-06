# Set variables
dir1="/home/inghirami/echo-qgp-mhd-test/outr0001/postproc/readx/"
dir2="/home/inghirami/echo-qgp-mhd-test-NOB/outr0001/postproc/readx/"
title1="With B"
title2="Without B"


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
var="vx"
iniframe=1
endframe=201

# Set the axes
set xlabel 'x (fm)'
#set ylabel "{/Symbol h}"
set ylabel 'y (fm)'
set cblabel "V_x (c units)"
set xrange [-15:15]
set yrange [-15:15]


# Set the color palette
# This section should be used only for quantities which decrease with time (en. density, temperature, etc)
#i0 = sprintf(dir1.var.'%04.0f.dat',iniframe)
#stats i0 using 3
#stmin=STATS_min
#stmax=STATS_max
#step=stmax-stmin
stmin=-1
stmax=1
step=stmax-stmin
set cbrange[stmin:stmax]
#set palette model RGB defined (stmin 'black', stmin+step/6 'blue', stmin+step/3 'cyan' , stmin+step/2 'green' , stmin+step*2/3. 'yellow', stmin+step*5/6. 'red', stmax 'purple')
set palette model RGB defined (stmin 'black', stmin+step/2 'red' , stmax 'yellow')

do for [i=iniframe:endframe] {
   keytitle="z=0 - t=".sprintf("%5.2f",(i-1)*0.02+0.5)." fm/c"
   infile1 = sprintf(dir1.var.'%04.0f.dat',i)
   infile2 = sprintf(dir2.var.'%04.0f.dat',i)
   outfile = sprintf(var.'%04.0f.png',i)
   set output outfile
   set multiplot layout 1,2 title keytitle
   set title title1
   splot infile1 u 1:2:3 notitle
   set title title2
   splot infile2 u 1:2:3 notitle
   unset multiplot
   print "Frame ".i." done"
}                

