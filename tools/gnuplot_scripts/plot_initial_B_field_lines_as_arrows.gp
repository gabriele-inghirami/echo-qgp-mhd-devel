#script date: 8/11/2018

set term png enh font "Helvetica, 16"

set xlabel "x (fm)"
set ylabel "y (fm)"

inputfile="Barrows12.dat"
plot_description_title="(Au+Au \\@ 200GeV, b=12 fm, {/Symbol h}=0, {/Symbol t}=0.4 fm, units: m^2_{{/Symbol p}_0})"

#to go from (GeV/fm^3)^1/2 to units of pion mass squared
conv_factor=4.81127
#conv_factor=1

iH1=5
iH2=8
iHC1=11
iHC2=14
iHt=17

stats inputfile using iH1
H1max=STATS_max
stats inputfile using iH2
H2max=STATS_max
stats inputfile using iHC1
HC1max=STATS_max
stats inputfile using iHC2
HC2max=STATS_max
stats inputfile using iHt
Htmax=STATS_max

set size square
set xrange [-15:15]
set yrange [-15:15]

set out "Bphi1_arrows_new.png"
set title "B_{/Symbol f} 1st ion ".plot_description_title
set palette defined (0 "white", H1max*conv_factor/4 "orange", H1max*conv_factor/2 "red", 3*H1max*conv_factor/4 "blue", H1max*conv_factor "green")
plot inputfile u (($1)-($3)/2.):(($2)-($4)/2.):3:4:($5)*conv_factor with vectors head size 0.5,20,60 filled lc palette lt 3 notitle

set out "Bphi2_arrows_new.png"
set title "B_{/Symbol f} 2nd ion ".plot_description_title
set palette defined (0 "white", H2max*conv_factor/4 "orange", H2max*conv_factor/2 "red", 3*H2max*conv_factor/4 "blue", H2max*conv_factor "green")
plot inputfile u (($1)-($6)/2.):(($2)-($7)/2.):6:7:($8)*conv_factor with vectors head size 0.5,20,60 filled lc palette lt 3 notitle

set out "Br1_arrows_new.png"
set title "B_r 1st ion ".plot_description_title
set palette defined (0 "white", HC1max*conv_factor/4 "orange", HC1max*conv_factor/2 "red", 3*HC1max*conv_factor/4 "blue", HC1max*conv_factor "green")
plot inputfile u (($1)-($9)/2.):(($2)-($10)/2.):9:10:($11)*conv_factor with vectors head size 0.5,20,60 filled lc palette lt 3 notitle

set out "Br2_arrows_new.png"
set title "B_r 2nd ion ".plot_description_title
set palette defined (0 "white", HC2max*conv_factor/4 "orange", HC2max*conv_factor/2 "red", 3*HC2max*conv_factor/4 "blue", HC2max*conv_factor "green")
plot inputfile u (($1)-($12)/2.):(($2)-($13)/2.):12:13:($14)*conv_factor with vectors head size 0.5,20,60 filled lc palette lt 3 notitle

set out "Btot_arrows_new.png"
set title "B_tot ".plot_description_title
set palette defined (0 "white", Htmax*conv_factor/4 "orange", Htmax*conv_factor/2 "red", 3*Htmax*conv_factor/4 "blue", Htmax*conv_factor "green")
plot inputfile u (($1)-($15)/2.):(($2)-($16)/2.):15:16:($17)*conv_factor with vectors head size 0.5,20,60 filled lc palette lt 3 notitle

