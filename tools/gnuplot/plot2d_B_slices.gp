# This script makes 2D plots of the initial magnetic fields
# used to initialize ECHO-QGP
set term png enh font "Helvetica, 14" size 800,600

set xlabel "x (fm)"
set ylabel "y (fm)"
set cblabel "units of m^2_{{/Symbol p}_0}"

inputfile="B12_1_1_4gp"

index_min=19
index_max=95

st=0.2

xside=15
yside=15
zside=11.5

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

stats inputfile using iH1x
H1xmin=STATS_min
H1xmax=STATS_max
stats inputfile using iH2x
H2xmin=STATS_min
H2xmax=STATS_max
stats inputfile using iHx
Hxmin=STATS_min
Hxmax=STATS_max
stats inputfile using iH1Cx
H1Cxmin=STATS_min
H1Cxmax=STATS_max
stats inputfile using iH2Cx
H2Cxmin=STATS_min
H2Cxmax=STATS_max
stats inputfile using iHCx
HCxmin=STATS_min
HCxmax=STATS_max
stats inputfile using iHxtot
Hxtotmin=STATS_min
Hxtotmax=STATS_max
stats inputfile using iH1y
H1ymin=STATS_min
H1ymax=STATS_max
stats inputfile using iH2y
H2ymin=STATS_min
H2ymax=STATS_max
stats inputfile using iHy
Hymin=STATS_min
Hymax=STATS_max
stats inputfile using iH1Cy
H1Cymin=STATS_min
H1Cymax=STATS_max
stats inputfile using iH2Cy
H2Cymin=STATS_min
H2Cymax=STATS_max
stats inputfile using iHCy
HCymin=STATS_min
HCymax=STATS_max
stats inputfile using iHytot
Hytotmin=STATS_min
Hytotmax=STATS_max
stats inputfile using iH1Cz
H1Czmin=STATS_min
H1Czmax=STATS_max
stats inputfile using iH2Cz
H2Czmin=STATS_min
H2Czmax=STATS_max
stats inputfile using iHCz
HCzmin=STATS_min
HCzmax=STATS_max

do for[i=index_min:index_max] {


pos=-zside+(i+0.5)*st

posstring=sprintf( "-eta%+-5.2f", pos )
postitlestring=sprintf( "-{/Symbol h}=%+-5.2f", pos )

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
 
set out "Bxs1".posstring.".png"
set title "B_x (1) classic ".postitlestring
STATS_min=H1xmin
STATS_max=H1xmax
set cbrange [STATS_min:STATS_max]
set palette model RGB defined ( STATS_min 'black', 0.67*STATS_min 'blue', 0.33*STATS_min 'cyan' ,0 'green' ,0 'yellow', 0.33*STATS_max 'red', 0.67*STATS_max 'purple', STATS_max 'white' )
splot inputfile index i using 1:2:4 notitle with lines lt 1 lw 0


set out "Bxs2".posstring.".png"
set title "B_x (2) classic ".postitlestring
STATS_min=H2xmin
STATS_max=H2xmax
set cbrange [STATS_min:STATS_max]
set palette model RGB defined ( STATS_min 'black', 0.67*STATS_min 'blue', 0.33*STATS_min 'cyan' ,0 'green' ,0 'yellow', 0.33*STATS_max 'red', 0.67*STATS_max 'purple', STATS_max 'white' )
splot inputfile index i using 1:2:iH2x notitle with lines lt 1 lw 0


set out "Bxst".posstring.".png"
set title "B_x classic".postitlestring
STATS_min=Hxmin
STATS_max=Hxmax
set cbrange [STATS_min:STATS_max]
set palette model RGB defined ( STATS_min 'black', 0.67*STATS_min 'blue', 0.33*STATS_min 'cyan' ,0 'green' ,0 'yellow', 0.33*STATS_max 'red', 0.67*STATS_max 'purple', STATS_max 'white' )
splot inputfile index i using 1:2:iHx notitle with lines lt 1 lw 0


set out "Bxc1".posstring.".png"
set title "B_x (1) chiral ".postitlestring
STATS_min=H1Cxmin
STATS_max=H1Cxmax
set cbrange [STATS_min:STATS_max]
set palette model RGB defined ( STATS_min 'black', 0.67*STATS_min 'blue', 0.33*STATS_min 'cyan' ,0 'green' ,0 'yellow', 0.33*STATS_max 'red', 0.67*STATS_max 'purple', STATS_max 'white' )
splot inputfile index i using 1:2:iH1Cx notitle with lines lt 1 lw 0


set out "Bxc2".posstring.".png"
set title "B_x (2) chiral ".postitlestring
STATS_min=H2Cxmin
STATS_max=H2Cxmax
set cbrange [STATS_min:STATS_max]
set palette model RGB defined ( STATS_min 'black', 0.67*STATS_min 'blue', 0.33*STATS_min 'cyan' ,0 'green' ,0 'yellow', 0.33*STATS_max 'red', 0.67*STATS_max 'purple', STATS_max 'white' )
splot inputfile index i using 1:2:iH2Cx notitle with lines lt 1 lw 0


set out "Bxct".posstring.".png"
set title "B_x chiral".postitlestring
STATS_min=HCxmin
STATS_max=HCxmax
set cbrange [STATS_min:STATS_max]
set palette model RGB defined ( STATS_min 'black', 0.67*STATS_min 'blue', 0.33*STATS_min 'cyan' ,0 'green' ,0 'yellow', 0.33*STATS_max 'red', 0.67*STATS_max 'purple', STATS_max 'white' )
splot inputfile index i using 1:2:iHCx notitle with lines lt 1 lw 0


set out "Bxt".posstring.".png"
set title "B_x total".postitlestring
STATS_min=Hxtotmin
STATS_max=Hxtotmax
set cbrange [STATS_min:STATS_max]
set palette model RGB defined ( STATS_min 'black', 0.67*STATS_min 'blue', 0.33*STATS_min 'cyan' ,0 'green' ,0 'yellow', 0.33*STATS_max 'red', 0.67*STATS_max 'purple', STATS_max 'white' )
splot inputfile index i using 1:2:iHxtot notitle with lines lt 1 lw 0


set out "Bys1".posstring.".png"
set title "B_y (1) classic ".postitlestring
STATS_min=H1ymin
STATS_max=H1ymax
set cbrange [STATS_min:STATS_max]
set palette model RGB defined ( STATS_min 'black', 0.67*STATS_min 'blue', 0.33*STATS_min 'cyan' ,0 'green' ,0 'yellow', 0.33*STATS_max 'red', 0.67*STATS_max 'purple', STATS_max 'white' )
splot inputfile index i using 1:2:iH1y notitle with lines lt 1 lw 0


set out "Bys2".posstring.".png"
set title "B_y (2) classic ".postitlestring
STATS_min=H2ymin
STATS_max=H2ymax
set cbrange [STATS_min:STATS_max]
set palette model RGB defined ( STATS_min 'black', 0.67*STATS_min 'blue', 0.33*STATS_min 'cyan' ,0 'green' ,0 'yellow', 0.33*STATS_max 'red', 0.67*STATS_max 'purple', STATS_max 'white' )
splot inputfile index i using 1:2:iH2y notitle with lines lt 1 lw 0


set out "Byst".posstring.".png"
set title "B_y classic".postitlestring
STATS_min=Hymin
STATS_max=Hymax
set cbrange [STATS_min:STATS_max]
set palette model RGB defined ( STATS_min 'black', 0.67*STATS_min 'blue', 0.33*STATS_min 'cyan' ,0 'green' ,0 'yellow', 0.33*STATS_max 'red', 0.67*STATS_max 'purple', STATS_max 'white' )
splot inputfile index i using 1:2:iHy notitle with lines lt 1 lw 0


set out "Byc1".posstring.".png"
set title "B_y (1) chiral ".postitlestring
STATS_min=H1Cymin
STATS_max=H1Cymax
set cbrange [STATS_min:STATS_max]
set palette model RGB defined ( STATS_min 'black', 0.67*STATS_min 'blue', 0.33*STATS_min 'cyan' ,0 'green' ,0 'yellow', 0.33*STATS_max 'red', 0.67*STATS_max 'purple', STATS_max 'white' )
splot inputfile index i using 1:2:iH1Cy notitle with lines lt 1 lw 0


set out "Byc2".posstring.".png"
set title "B_y (2) chiral ".postitlestring
STATS_min=H2Cymin
STATS_max=H2Cymax
set cbrange [STATS_min:STATS_max]
set palette model RGB defined ( STATS_min 'black', 0.67*STATS_min 'blue', 0.33*STATS_min 'cyan' ,0 'green' ,0 'yellow', 0.33*STATS_max 'red', 0.67*STATS_max 'purple', STATS_max 'white' )
splot inputfile index i using 1:2:iH2Cy notitle with lines lt 1 lw 0


set out "Byct".posstring.".png"
set title "B_y chiral".postitlestring
STATS_min=HCymin
STATS_max=HCymax
set cbrange [STATS_min:STATS_max]
set palette model RGB defined ( STATS_min 'black', 0.67*STATS_min 'blue', 0.33*STATS_min 'cyan' ,0 'green' ,0 'yellow', 0.33*STATS_max 'red', 0.67*STATS_max 'purple', STATS_max 'white' )
splot inputfile index i using 1:2:iHCy notitle with lines lt 1 lw 0


set out "Byt".posstring.".png"
set title "B_y total".postitlestring
STATS_min=Hytotmin
STATS_max=Hytotmax
set cbrange [STATS_min:STATS_max]
set palette model RGB defined ( STATS_min 'black', 0.67*STATS_min 'blue', 0.33*STATS_min 'cyan' ,0 'green' ,0 'yellow', 0.33*STATS_max 'red', 0.67*STATS_max 'purple', STATS_max 'white' )
splot inputfile index i using 1:2:iHytot notitle with lines lt 1 lw 0


set out "Bzc1".posstring.".png"
set title "B_z (1) chiral ".postitlestring
STATS_min=H1Czmin
STATS_max=H1Czmax
set cbrange [STATS_min:STATS_max]
set palette model RGB defined ( STATS_min 'black', 0.67*STATS_min 'blue', 0.33*STATS_min 'cyan' ,0 'green' ,0 'yellow', 0.33*STATS_max 'red', 0.67*STATS_max 'purple', STATS_max 'white' )
splot inputfile index i using 1:2:iH1Cz notitle with lines lt 1 lw 0


set out "Bzc2".posstring.".png"
set title "B_z (2) chiral ".postitlestring
STATS_min=H2Czmin
STATS_max=H2Czmax
set cbrange [STATS_min:STATS_max]
set palette model RGB defined ( STATS_min 'black', 0.67*STATS_min 'blue', 0.33*STATS_min 'cyan' ,0 'green' ,0 'yellow', 0.33*STATS_max 'red', 0.67*STATS_max 'purple', STATS_max 'white' )
splot inputfile index i using 1:2:iH2Cz notitle with lines lt 1 lw 0


set out "Bzt".posstring.".png"
set title "B_z chiral ".postitlestring
STATS_min=HCzmin
STATS_max=HCzmax
set cbrange [STATS_min:STATS_max]
set palette model RGB defined ( STATS_min 'black', 0.67*STATS_min 'blue', 0.33*STATS_min 'cyan' ,0 'green' ,0 'yellow', 0.33*STATS_max 'red', 0.67*STATS_max 'purple', STATS_max 'white' )
splot inputfile index i using 1:2:iHCz notitle with lines lt 1 lw 0


}
