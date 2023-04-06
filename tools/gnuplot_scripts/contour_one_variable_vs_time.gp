# Set variables
echo_qgp_dir="../../"
dir1=echo_qgp_dir."outr0001/postproc/readx/"

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
var="B_x"
iniframe=1
endframe=201

# Set the axes
set xlabel 'x (fm)'
#set ylabel "{/Symbol h}"
set ylabel 'y (fm)'
set cblabel "B_x (c units)"
set xrange [-15:15]
set yrange [-15:15]


# Set the color palette
#i0 = sprintf(dir1.var.'%04.0f.dat',iniframe)
#stats i0 using 3
#stmin=STATS_min
#stmax=STATS_max
stmin=-0.1
stmax=0.1
step=stmax-stmin
set cbrange[stmin:stmax]
#set palette model RGB defined (stmin 'black', stmin+step/6 'blue', stmin+step/3 'cyan' , stmin+step/2 'green' , stmin+step*2/3. 'yellow', stmin+step*5/6. 'red', stmax 'purple')
set palette model RGB defined (stmin 'black', stmin+step/2 'red' , stmax 'yellow')

do for [i=iniframe:endframe] {
   keytitle="z=0 - t=".sprintf("%5.2f",(i-1)*0.02+0.5)." fm/c"
   infile1 = sprintf(dir1.var.'%04.0f.dat',i)
   outfile = sprintf(var.'%04.0f.png',i)
   set output outfile
   set title keytitle
   splot infile1 u 1:2:3 notitle
   print "Frame ".i." done"
}                

