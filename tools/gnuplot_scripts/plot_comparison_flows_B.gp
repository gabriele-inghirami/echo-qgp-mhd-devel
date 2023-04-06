# sample script to compare v1(y), v2(y), v1(pT), v2(pT)
# of positive pions for three different simulations
set term pos eps color enh font "Helvetic, 20"

echo_qgp_dir="../../"
d1=echo_qgp_dir."outr0009_run1/spectra/"
d2=echo_qgp_dir."outr0009_run2/spectra/"
d3=echo_qgp_dir."outr0009_run3/spectra/"

t1="run 1"
t2="run 2"
t3="run 3"

filerap="rap_spectra__00000211_+_Pion_plus.txt"
filept="spectra_v2__00000211_+_Pion_plus.txt"


set xlabel "y (rapidity)"

set out "V1_vs_Y.eps"
set ylabel "v_1"
set key top left
plot d1.filerap u 1:5 w l lw 2 t t1, d2.filerap u 1:5 w l lw 2 t t2, d3.filerap u 1:5 w l lw 2 t t3

set out "V2_vs_Y.eps"
set key top right
set ylabel "v_2"
plot d1.filerap u 1:6 w l lw 2 t t1, d2.filerap u 1:6 w l lw 2 t t2, d3.filerap u 1:6 w l lw 2 t t3


set xlabel "pT [GeV]"

set out "V1_vs_pT.eps"
set ylabel "v_1"
set key top right
plot d1.filept i 21 u 1:4 w l lw 2 t t1, d2.filept i 21 u 1:4 w l lw 2 t t2, d3.filept i 21 u 1:4 w l lw 2 t t3

set out "V2_vs_pT.eps"
set key top left
set ylabel "v_2"
plot d1.filept i 21 u 1:5 w l lw 2 t t1, d2.filept i 21 u 1:5 w l lw 2 t t2, d3.filept i 21 u 1:5 w l lw 2 t t3
