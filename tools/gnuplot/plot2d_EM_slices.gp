# This script makes 2D plots of the initial magnetic and electric
# fields used to initialize ECHO-QGP.
# Each plot uses its own maximum and minimum values
set term png enh font "Helvetica, 14" size 800,600

set xlabel "x (fm)"
set ylabel "y (fm)"
set cblabel "units of m^2_{{/Symbol p}_0}"

inputfile="EM_zplu1_4gp"

index_min=0
index_max=0

st=0.0 #the dz step

xside=20
yside=20
zside=0. #the left edge of the z axis 

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

stats inputfile using iE1x
E1xmin=STATS_min
E1xmax=STATS_max
stats inputfile using iE2x
E2xmin=STATS_min
E2xmax=STATS_max
stats inputfile using iEx
Exmin=STATS_min
Exmax=STATS_max
stats inputfile using iE1Cx
E1Cxmin=STATS_min
E1Cxmax=STATS_max
stats inputfile using iE2Cx
E2Cxmin=STATS_min
E2Cxmax=STATS_max
stats inputfile using iECx
ECxmin=STATS_min
ECxmax=STATS_max
stats inputfile using iExtot
Extotmin=STATS_min
Extotmax=STATS_max
stats inputfile using iE1y
E1ymin=STATS_min
E1ymax=STATS_max
stats inputfile using iE2y
E2ymin=STATS_min
E2ymax=STATS_max
stats inputfile using iEy
Eymin=STATS_min
Eymax=STATS_max
stats inputfile using iE1Cy
E1Cymin=STATS_min
E1Cymax=STATS_max
stats inputfile using iE2Cy
E2Cymin=STATS_min
E2Cymax=STATS_max
stats inputfile using iECy
ECymin=STATS_min
ECymax=STATS_max
stats inputfile using iEytot
Eytotmin=STATS_min
Eytotmax=STATS_max
stats inputfile using iE1z
E1zmin=STATS_min
E1zmax=STATS_max
stats inputfile using iE2z
E2zmin=STATS_min
E2zmax=STATS_max
stats inputfile using iEz
Ezmin=STATS_min
Ezmax=STATS_max

do for[i=index_min:index_max] {


#pos=-zside+(i+0.5)*st
pos=+1

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

if(H1xmin!=H1xmax) { set out "Bxs1".posstring.".png"; set title "B_x (1) classic ".postitlestring; STATS_min=H1xmin; STATS_max=H1xmax; set cbrange [STATS_min:STATS_max]; set palette model RGB defined ( STATS_min 'black', 0.67*STATS_min 'blue', 0.33*STATS_min 'cyan' ,0 'green' ,0 'yellow', 0.33*STATS_max 'red', 0.67*STATS_max 'purple', STATS_max 'white' ); splot inputfile index i using 1:2:4 notitle with lines lt 1 lw 0; }

if(H2xmin!=H2xmax) { set out "Bxs2".posstring.".png"; set title "B_x (2) classic ".postitlestring; STATS_min=H2xmin; STATS_max=H2xmax; set cbrange [STATS_min:STATS_max]; set palette model RGB defined ( STATS_min 'black', 0.67*STATS_min 'blue', 0.33*STATS_min 'cyan' ,0 'green' ,0 'yellow', 0.33*STATS_max 'red', 0.67*STATS_max 'purple', STATS_max 'white' ); splot inputfile index i using 1:2:iH2x notitle with lines lt 1 lw 0; }


if(Hxmin!=Hxmax) { set out "Bxst".posstring.".png"; set title "B_x classic".postitlestring; STATS_min=Hxmin; STATS_max=Hxmax; set cbrange [STATS_min:STATS_max]; set palette model RGB defined ( STATS_min 'black', 0.67*STATS_min 'blue', 0.33*STATS_min 'cyan' ,0 'green' ,0 'yellow', 0.33*STATS_max 'red', 0.67*STATS_max 'purple', STATS_max 'white' ); splot inputfile index i using 1:2:iHx notitle with lines lt 1 lw 0; }


if(H1Cxmin!=H1Cxmax) { set out "Bxc1".posstring.".png"; set title "B_x (1) chiral ".postitlestring; STATS_min=H1Cxmin; STATS_max=H1Cxmax; set cbrange [STATS_min:STATS_max]; set palette model RGB defined ( STATS_min 'black', 0.67*STATS_min 'blue', 0.33*STATS_min 'cyan' ,0 'green' ,0 'yellow', 0.33*STATS_max 'red', 0.67*STATS_max 'purple', STATS_max 'white' ); splot inputfile index i using 1:2:iH1Cx notitle with lines lt 1 lw 0; }


if(H2Cxmin!=H2Cxmax) { set out "Bxc2".posstring.".png"; set title "B_x (2) chiral ".postitlestring; STATS_min=H2Cxmin; STATS_max=H2Cxmax; set cbrange [STATS_min:STATS_max]; set palette model RGB defined ( STATS_min 'black', 0.67*STATS_min 'blue', 0.33*STATS_min 'cyan' ,0 'green' ,0 'yellow', 0.33*STATS_max 'red', 0.67*STATS_max 'purple', STATS_max 'white' ); splot inputfile index i using 1:2:iH2Cx notitle with lines lt 1 lw 0; }


if(HCxmin!=HCxmax) { set out "Bxct".posstring.".png"; set title "B_x chiral".postitlestring; STATS_min=HCxmin; STATS_max=HCxmax; set cbrange [STATS_min:STATS_max]; set palette model RGB defined ( STATS_min 'black', 0.67*STATS_min 'blue', 0.33*STATS_min 'cyan' ,0 'green' ,0 'yellow', 0.33*STATS_max 'red', 0.67*STATS_max 'purple', STATS_max 'white' ); splot inputfile index i using 1:2:iHCx notitle with lines lt 1 lw 0; }


if(Hxtotmin!=Hxtotmax) { set out "Bxt".posstring.".png"; set title "B_x total".postitlestring; STATS_min=Hxtotmin; STATS_max=Hxtotmax; set cbrange [STATS_min:STATS_max]; set palette model RGB defined ( STATS_min 'black', 0.67*STATS_min 'blue', 0.33*STATS_min 'cyan' ,0 'green' ,0 'yellow', 0.33*STATS_max 'red', 0.67*STATS_max 'purple', STATS_max 'white' ); splot inputfile index i using 1:2:iHxtot notitle with lines lt 1 lw 0; }


if(H1ymin!=H1ymax) { set out "Bys1".posstring.".png"; set title "B_y (1) classic ".postitlestring; STATS_min=H1ymin; STATS_max=H1ymax; set cbrange [STATS_min:STATS_max]; set palette model RGB defined ( STATS_min 'black', 0.67*STATS_min 'blue', 0.33*STATS_min 'cyan' ,0 'green' ,0 'yellow', 0.33*STATS_max 'red', 0.67*STATS_max 'purple', STATS_max 'white' ); splot inputfile index i using 1:2:iH1y notitle with lines lt 1 lw 0; }


if(H2ymin!=H2ymax) { set out "Bys2".posstring.".png"; set title "B_y (2) classic ".postitlestring; STATS_min=H2ymin; STATS_max=H2ymax; set cbrange [STATS_min:STATS_max]; set palette model RGB defined ( STATS_min 'black', 0.67*STATS_min 'blue', 0.33*STATS_min 'cyan' ,0 'green' ,0 'yellow', 0.33*STATS_max 'red', 0.67*STATS_max 'purple', STATS_max 'white' ); splot inputfile index i using 1:2:iH2y notitle with lines lt 1 lw 0; }


if(Hymin!=Hymax) { set out "Byst".posstring.".png"; set title "B_y classic".postitlestring; STATS_min=Hymin; STATS_max=Hymax; set cbrange [STATS_min:STATS_max]; set palette model RGB defined ( STATS_min 'black', 0.67*STATS_min 'blue', 0.33*STATS_min 'cyan' ,0 'green' ,0 'yellow', 0.33*STATS_max 'red', 0.67*STATS_max 'purple', STATS_max 'white' ); splot inputfile index i using 1:2:iHy notitle with lines lt 1 lw 0; }


if(H1Cymin!=H1Cymax) { set out "Byc1".posstring.".png"; set title "B_y (1) chiral ".postitlestring; STATS_min=H1Cymin; STATS_max=H1Cymax; set cbrange [STATS_min:STATS_max]; set palette model RGB defined ( STATS_min 'black', 0.67*STATS_min 'blue', 0.33*STATS_min 'cyan' ,0 'green' ,0 'yellow', 0.33*STATS_max 'red', 0.67*STATS_max 'purple', STATS_max 'white' ); splot inputfile index i using 1:2:iH1Cy notitle with lines lt 1 lw 0; }


if(H2Cymin!=H2Cymax) { set out "Byc2".posstring.".png"; set title "B_y (2) chiral ".postitlestring; STATS_min=H2Cymin; STATS_max=H2Cymax; set cbrange [STATS_min:STATS_max]; set palette model RGB defined ( STATS_min 'black', 0.67*STATS_min 'blue', 0.33*STATS_min 'cyan' ,0 'green' ,0 'yellow', 0.33*STATS_max 'red', 0.67*STATS_max 'purple', STATS_max 'white' ); splot inputfile index i using 1:2:iH2Cy notitle with lines lt 1 lw 0; }


if(HCymin!=HCymax) { set out "Byct".posstring.".png"; set title "B_y chiral".postitlestring; STATS_min=HCymin; STATS_max=HCymax; set cbrange [STATS_min:STATS_max]; set palette model RGB defined ( STATS_min 'black', 0.67*STATS_min 'blue', 0.33*STATS_min 'cyan' ,0 'green' ,0 'yellow', 0.33*STATS_max 'red', 0.67*STATS_max 'purple', STATS_max 'white' ); splot inputfile index i using 1:2:iHCy notitle with lines lt 1 lw 0; }


if(Hytotmin!=Hytotmax) { set out "Byt".posstring.".png"; set title "B_y total".postitlestring; STATS_min=Hytotmin; STATS_max=Hytotmax; set cbrange [STATS_min:STATS_max]; set palette model RGB defined ( STATS_min 'black', 0.67*STATS_min 'blue', 0.33*STATS_min 'cyan' ,0 'green' ,0 'yellow', 0.33*STATS_max 'red', 0.67*STATS_max 'purple', STATS_max 'white' ); splot inputfile index i using 1:2:iHytot notitle with lines lt 1 lw 0; }


if(H1Czmin!=H1Czmax) { set out "Bzc1".posstring.".png"; set title "B_z (1) chiral ".postitlestring; STATS_min=H1Czmin; STATS_max=H1Czmax; set cbrange [STATS_min:STATS_max]; set palette model RGB defined ( STATS_min 'black', 0.67*STATS_min 'blue', 0.33*STATS_min 'cyan' ,0 'green' ,0 'yellow', 0.33*STATS_max 'red', 0.67*STATS_max 'purple', STATS_max 'white' ); splot inputfile index i using 1:2:iH1Cz notitle with lines lt 1 lw 0; }


if(H2Czmin!=H2Czmax) { set out "Bzc2".posstring.".png"; set title "B_z (2) chiral ".postitlestring; STATS_min=H2Czmin; STATS_max=H2Czmax; set cbrange [STATS_min:STATS_max]; set palette model RGB defined ( STATS_min 'black', 0.67*STATS_min 'blue', 0.33*STATS_min 'cyan' ,0 'green' ,0 'yellow', 0.33*STATS_max 'red', 0.67*STATS_max 'purple', STATS_max 'white' ); splot inputfile index i using 1:2:iH2Cz notitle with lines lt 1 lw 0; }


if(HCzmin!=HCzmax) { set out "Bzt".posstring.".png"; set title "B_z chiral ".postitlestring; STATS_min=HCzmin; STATS_max=HCzmax; set cbrange [STATS_min:STATS_max]; set palette model RGB defined ( STATS_min 'black', 0.67*STATS_min 'blue', 0.33*STATS_min 'cyan' ,0 'green' ,0 'yellow', 0.33*STATS_max 'red', 0.67*STATS_max 'purple', STATS_max 'white' ); splot inputfile index i using 1:2:iHCz notitle with lines lt 1 lw 0; }


if(E1xmin!=E1xmax) { set out "Exs1".posstring.".png"; set title "E_x (1) classic ".postitlestring; STATS_min=E1xmin; STATS_max=E1xmax; set cbrange [STATS_min:STATS_max]; set palette model RGB defined ( STATS_min 'black', 0.67*STATS_min 'blue', 0.33*STATS_min 'cyan' ,0 'green' ,0 'yellow', 0.33*STATS_max 'red', 0.67*STATS_max 'purple', STATS_max 'white' ); splot inputfile index i using 1:2:iE1x notitle with lines lt 1 lw 0; }


if(E2xmin!=E2xmax) { set out "Exs2".posstring.".png"; set title "E_x (2) classic ".postitlestring; STATS_min=E2xmin; STATS_max=E2xmax; set cbrange [STATS_min:STATS_max]; set palette model RGB defined ( STATS_min 'black', 0.67*STATS_min 'blue', 0.33*STATS_min 'cyan' ,0 'green' ,0 'yellow', 0.33*STATS_max 'red', 0.67*STATS_max 'purple', STATS_max 'white' ); splot inputfile index i using 1:2:iE2x notitle with lines lt 1 lw 0; }


if(Exmin!=Exmax) { set out "Exst".posstring.".png"; set title "E_x classic".postitlestring; STATS_min=Exmin; STATS_max=Exmax; set cbrange [STATS_min:STATS_max]; set palette model RGB defined ( STATS_min 'black', 0.67*STATS_min 'blue', 0.33*STATS_min 'cyan' ,0 'green' ,0 'yellow', 0.33*STATS_max 'red', 0.67*STATS_max 'purple', STATS_max 'white' ); splot inputfile index i using 1:2:iEx notitle with lines lt 1 lw 0; }


if(E1Cxmin!=E1Cxmax) { set out "Exc1".posstring.".png"; set title "E_x (1) chiral ".postitlestring; STATS_min=E1Cxmin; STATS_max=E1Cxmax; set cbrange [STATS_min:STATS_max]; set palette model RGB defined ( STATS_min 'black', 0.67*STATS_min 'blue', 0.33*STATS_min 'cyan' ,0 'green' ,0 'yellow', 0.33*STATS_max 'red', 0.67*STATS_max 'purple', STATS_max 'white' ); splot inputfile index i using 1:2:iE1Cx notitle with lines lt 1 lw 0; }


if(E2Cxmin!=E2Cxmax) { set out "Exc2".posstring.".png"; set title "E_x (2) chiral ".postitlestring; STATS_min=E2Cxmin; STATS_max=E2Cxmax; set cbrange [STATS_min:STATS_max]; set palette model RGB defined ( STATS_min 'black', 0.67*STATS_min 'blue', 0.33*STATS_min 'cyan' ,0 'green' ,0 'yellow', 0.33*STATS_max 'red', 0.67*STATS_max 'purple', STATS_max 'white' ); splot inputfile index i using 1:2:iE2Cx notitle with lines lt 1 lw 0; }


if(ECxmin!=ECxmax) { set out "Exct".posstring.".png"; set title "E_x chiral".postitlestring; STATS_min=ECxmin; STATS_max=ECxmax; set cbrange [STATS_min:STATS_max]; set palette model RGB defined ( STATS_min 'black', 0.67*STATS_min 'blue', 0.33*STATS_min 'cyan' ,0 'green' ,0 'yellow', 0.33*STATS_max 'red', 0.67*STATS_max 'purple', STATS_max 'white' ); splot inputfile index i using 1:2:iECx notitle with lines lt 1 lw 0; }


if(Extotmin!=Extotmax) { set out "Ext".posstring.".png"; set title "E_x total".postitlestring; STATS_min=Extotmin; STATS_max=Extotmax; set cbrange [STATS_min:STATS_max]; set palette model RGB defined ( STATS_min 'black', 0.67*STATS_min 'blue', 0.33*STATS_min 'cyan' ,0 'green' ,0 'yellow', 0.33*STATS_max 'red', 0.67*STATS_max 'purple', STATS_max 'white' ); splot inputfile index i using 1:2:iExtot notitle with lines lt 1 lw 0; }


if(E1ymin!=E1ymax) { set out "Eys1".posstring.".png"; set title "E_y (1) classic ".postitlestring; STATS_min=E1ymin; STATS_max=E1ymax; set cbrange [STATS_min:STATS_max]; set palette model RGB defined ( STATS_min 'black', 0.67*STATS_min 'blue', 0.33*STATS_min 'cyan' ,0 'green' ,0 'yellow', 0.33*STATS_max 'red', 0.67*STATS_max 'purple', STATS_max 'white' ); splot inputfile index i using 1:2:iE1y notitle with lines lt 1 lw 0; }


if(E2ymin!=E2ymax) { set out "Eys2".posstring.".png"; set title "E_y (2) classic ".postitlestring; STATS_min=E2ymin; STATS_max=E2ymax; set cbrange [STATS_min:STATS_max]; set palette model RGB defined ( STATS_min 'black', 0.67*STATS_min 'blue', 0.33*STATS_min 'cyan' ,0 'green' ,0 'yellow', 0.33*STATS_max 'red', 0.67*STATS_max 'purple', STATS_max 'white' ); splot inputfile index i using 1:2:iE2y notitle with lines lt 1 lw 0; }


if(Eymin!=Eymax) { set out "Eyst".posstring.".png"; set title "E_y classic".postitlestring; STATS_min=Eymin; STATS_max=Eymax; set cbrange [STATS_min:STATS_max]; set palette model RGB defined ( STATS_min 'black', 0.67*STATS_min 'blue', 0.33*STATS_min 'cyan' ,0 'green' ,0 'yellow', 0.33*STATS_max 'red', 0.67*STATS_max 'purple', STATS_max 'white' ); splot inputfile index i using 1:2:iEy notitle with lines lt 1 lw 0; }


if(E1Cymin!=E1Cymax) { set out "Eyc1".posstring.".png"; set title "E_y (1) chiral ".postitlestring; STATS_min=E1Cymin; STATS_max=E1Cymax; set cbrange [STATS_min:STATS_max]; set palette model RGB defined ( STATS_min 'black', 0.67*STATS_min 'blue', 0.33*STATS_min 'cyan' ,0 'green' ,0 'yellow', 0.33*STATS_max 'red', 0.67*STATS_max 'purple', STATS_max 'white' ); splot inputfile index i using 1:2:iE1Cy notitle with lines lt 1 lw 0; }


if(E2Cymin!=E2Cymax) { set out "Eyc2".posstring.".png"; set title "E_y (2) chiral ".postitlestring; STATS_min=E2Cymin; STATS_max=E2Cymax; set cbrange [STATS_min:STATS_max]; set palette model RGB defined ( STATS_min 'black', 0.67*STATS_min 'blue', 0.33*STATS_min 'cyan' ,0 'green' ,0 'yellow', 0.33*STATS_max 'red', 0.67*STATS_max 'purple', STATS_max 'white' ); splot inputfile index i using 1:2:iE2Cy notitle with lines lt 1 lw 0; }


if(ECymin!=ECymax) { set out "Eyct".posstring.".png"; set title "E_y chiral".postitlestring; STATS_min=ECymin; STATS_max=ECymax; set cbrange [STATS_min:STATS_max]; set palette model RGB defined ( STATS_min 'black', 0.67*STATS_min 'blue', 0.33*STATS_min 'cyan' ,0 'green' ,0 'yellow', 0.33*STATS_max 'red', 0.67*STATS_max 'purple', STATS_max 'white' ); splot inputfile index i using 1:2:iECy notitle with lines lt 1 lw 0; }


if(Eytotmin!=Eytotmax) { set out "Eyt".posstring.".png"; set title "E_y total".postitlestring; STATS_min=Eytotmin; STATS_max=Eytotmax; set cbrange [STATS_min:STATS_max]; set palette model RGB defined ( STATS_min 'black', 0.67*STATS_min 'blue', 0.33*STATS_min 'cyan' ,0 'green' ,0 'yellow', 0.33*STATS_max 'red', 0.67*STATS_max 'purple', STATS_max 'white' ); splot inputfile index i using 1:2:iEytot notitle with lines lt 1 lw 0; }


if(E1zmin!=E1zmax) { set out "Ez1".posstring.".png"; set title "E_z (1) classic ".postitlestring; STATS_min=E1zmin; STATS_max=E1zmax; set cbrange [STATS_min:STATS_max]; if (STATS_min*STATS_max < 0){ set palette model RGB defined ( STATS_min 'black', 0.67*STATS_min 'blue', 0.33*STATS_min 'cyan' ,0 'green' ,0 'yellow', 0.33*STATS_max 'red', 0.67*STATS_max 'purple', STATS_max 'white' );}; else { set palette rgbformulae 33,13,10; }; splot inputfile index i using 1:2:iE1z notitle with lines lt 1 lw 0; }


if(E2zmin!=E2zmax) { set out "Ez2".posstring.".png"; set title "E_z (2) classic ".postitlestring; STATS_min=E2zmin; STATS_max=E2zmax; set cbrange [STATS_min:STATS_max]; if (STATS_min*STATS_max < 0) { set palette model RGB defined ( STATS_min 'black', 0.67*STATS_min 'blue', 0.33*STATS_min 'cyan' ,0 'green' ,0 'yellow', 0.33*STATS_max 'red', 0.67*STATS_max 'purple', STATS_max 'white' ); }; else { set palette rgbformulae 33,13,10; }; splot inputfile index i using 1:2:iE2z notitle with lines lt 1 lw 0; }


if(Ezmin!=Ezmax) { set out "Ezt".posstring.".png"; set title "E_z classic ".postitlestring; STATS_min=Ezmin; STATS_max=Ezmax; set cbrange [STATS_min:STATS_max]; if (STATS_min*STATS_max <=0) { set palette model RGB defined ( STATS_min 'black', 0.67*STATS_min 'blue', 0.33*STATS_min 'cyan' ,0 'green' ,0 'yellow', 0.33*STATS_max 'red', 0.67*STATS_max 'purple', STATS_max 'white' );} else { set palette rgbformulae 33,13,10; }; splot inputfile index i using 1:2:iEz notitle with lines lt 1 lw 0; }

}
