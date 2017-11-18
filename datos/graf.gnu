#!/usr/bin/gnuplot

reset

set title
set grid x,y
set samples 20000

#################################################################
#
# F(x)=exp(a)*exp(b*x)*exp(c*x**3) -> T2,  pomocna f(x)=log(F(x))
# b=-2/T_2
#
# G(x)=d*abs(1-exp(-e*x)) -> T1
# e=1/T_1
#
# rT=1/T,  rT_lin=alfalin*x linearni, rT=alfa*x+beta afinni
# u T_1 sedi rT_1lin
# u T_2 sedi lip rT_2
#
#
# fitlog:
#
# alfalin_1       = 0.00102649       +/- 1.406e-05    (1.369%)
#
# alfa_2          = 0.00132344       +/- 7.168e-05    (5.416%)  ms-1 mM-1
# beta_2          = -0.00671824      +/- 0.003804     (56.62%)  ms-1
#
# a_1             = -3.03358         +/- 0.005102     (0.1682%)
# b_1             = -0.0324671       +/- 0.0009456    (2.913%)
# c_1             = -5.52841e-05     +/- 2.23e-06     (4.034%)
#
# a_2             = -3.07068         +/- 0.008059     (0.2624%)
# b_2             = -0.0681121       +/- 0.00176      (2.583%)
# c_2             = -5.16632e-05     +/- 5.653e-06    (10.94%)
#
# a_3             = -2.97748         +/- 0.008069     (0.271%)
# b_3             = -0.115878        +/- 0.001983     (1.711%)
# c_3             = -1.89576e-05     +/- 8.079e-06    (42.62%)
#
# a_4             = -2.9715          +/- 0.008648     (0.291%)
# b_4             = -0.146674        +/- 0.002424     (1.652%)
# c_4             = -6.24124e-06     +/- 1.267e-05    (203%)
#
# a_5             = -3.34766         +/- 0.0123       (0.3673%)
# b_5             = -0.204936        +/- 0.004514     (2.203%)
# c_5             = 5.41352e-05      +/- 3.805e-05    (70.29%)
#
#
# d_1             = 0.0380539        +/- 0.0001326    (0.3485%)
# e_1             = 0.0157552        +/- 6.129e-05    (0.389%)
#
# d_2             = 0.0342925        +/- 0.0003924    (1.144%)
# e_2             = 0.0314669        +/- 0.0003828    (1.216%)
#
# d_3             = 0.0383184        +/- 0.0002809    (0.7329%)
# e_3             = 0.0488528        +/- 0.000376     (0.7696%)
#
# d_4             = 0.0388043        +/- 0.0004486    (1.156%)
# e_4             = 0.063935         +/- 0.0008222    (1.286%)
#
# d_5             = 0.0278024        +/- 0.0002637    (0.9484%)
# e_5             = 0.0844645        +/- 0.0009266    (1.097%)
#
#################################################################

########
## T2 ##
########

f_1(x)=a_1+b_1*x+c_1*x**3
F_1(x)=exp(a_1)*exp(b_1*x+c_1*x**3)
fit f_1(x) 'T2_v1.txt' u 1:(log($2)) via a_1, b_1, c_1
fit F_1(x) 'T2_v1.txt' via a_1, b_1, c_1

f_2(x)=a_2+b_2*x+c_2*x**3
F_2(x)=exp(a_2)*exp(b_2*x+c_2*x**3)
fit f_2(x) 'T2_v2.txt' u 1:(log($2)) via a_2, b_2, c_2
fit F_2(x) 'T2_v2.txt' via a_2, b_2, c_2

f_3(x)=a_3+b_3*x+c_3*x**3
F_3(x)=exp(a_3)*exp(b_3*x+c_3*x**3)
fit f_3(x) 'T2_v3.txt' u 1:(log($2)) via a_3, b_3, c_3
fit F_3(x) 'T2_v3.txt' via a_3, b_3, c_3

f_4(x)=a_4+b_4*x+c_4*x**3
F_4(x)=exp(a_4)*exp(b_4*x+c_4*x**3)
fit f_4(x) 'T2_v4.txt' u 1:(log($2)) via a_4, b_4, c_4
fit F_4(x) 'T2_v4.txt' via a_4, b_4, c_4

f_5(x)=a_5+b_5*x+c_5*x**3
F_5(x)=exp(a_5)*exp(b_5*x+c_5*x**3)
fit f_5(x) 'T2_v5.txt' u 1:(log($2)) via a_5, b_5, c_5
fit F_5(x) 'T2_v5.txt' via a_5, b_5, c_5

rT_1lin(x)=alfalin_1*x
rT_1(x)=alfa_1*x+beta_1
fit rT_1lin(x) 'koncentrace' u 2:7 via alfalin_1
fit rT_1(x) 'koncentrace' u 2:7 via alfa_1, beta_1


########
## T1 ##
########

G_1(x)=d_1*abs(1-2*exp(-e_1*x))
e_1=log(2)/20
fit G_1(x) 'T1_v1.txt' via d_1, e_1

G_2(x)=d_2*abs(1-2*exp(-e_2*x))
e_2=log(2)/20
fit G_2(x) 'T1_v2.txt' via d_2, e_2

G_3(x)=d_3*abs(1-2*exp(-e_3*x))
e_3=log(2)/20
fit G_3(x) 'T1_v3.txt' via d_3, e_3

G_4(x)=d_4*abs(1-2*exp(-e_4*x))
e_4=log(2)/20
fit G_4(x) 'T1_v4.txt' via d_4, e_4

G_5(x)=d_5*abs(1-2*exp(-e_5*x))
e_5=log(2)/20
fit G_5(x) 'T1_v5.txt' via d_5, e_5

rT_2lin(x)=alfalin_2*x
rT_2(x)=alfa_2*x+beta_2
fit rT_2lin(x) 'koncentrace' u 2:9 via alfalin_2
fit rT_2(x) 'koncentrace' u 2:9 via alfa_2, beta_2

##########
## plot ##
##########

set terminal epslatex size 18cm,12cm

set xlabel '$c$ (\si{\milli M})'
set xrange[0:90]

set output 'rT.tex'
set key top left
set yrange[0:0.11]
set ylabel '(\si{\per\ms})'
plot 'koncentrace' u 2:7:8 ps 3 t '$1/T_1$' w yerrorbars, rT_1lin(x) ls 1 lw 2 t 'lineární fit $1/T_1$' , '' u 2:9:($10) ls 2 ps 3 t '$1/T_2$' w yerrorbars, rT_2(x) ls 2 lw 2 t 'afinní fit $1/T_2$'

set output  'T.tex'
set key top right
set yrange[0:70]
set ylabel '(\si{\ms})'
plot 'koncentrace' u 2:3:4 ps 3 t '$T_1$' w yerrorbars, 1/rT_1lin(x) ls 1 lw 2 t 'lineární fit $1/T_1$' ,'' u 2:5:6 ls 2 ps 3 t '$T_2$' w yerrorbars, 1/rT_2(x) ls 2 lw 2 t 'afinní fit $1/T_2$'

set terminal epslatex size 18cm,24cm
set output 'g.tex'
set ylabel
set xlabel
unset ytics
set xtics
set multiplot layout 5, 2

set key bottom right
set yrange [0:0.04]
set xrange [0:270]
plot 'T1_v1.txt' ps 3 t 'IR, č. 1', G_1(x) ls 1 lw 2 notitle
set yrange [0:0.035]
set xrange [0:130]
plot 'T1_v2.txt' ps 3 t 'IR, č. 2', G_2(x) ls 1 lw 2 notitle
set yrange [0:0.04]
set xrange [0:110]
plot 'T1_v3.txt' ps 3 t 'IR, č. 3', G_3(x) ls 1 lw 2 notitle
set yrange [0:0.04]
set xrange [0:110]
plot 'T1_v4.txt' ps 3 t 'IR, č. 4', G_4(x) ls 1 lw 2 notitle
set yrange [0:0.03]
set xrange [0:90]
plot 'T1_v5.txt' ps 3 t 'IR, č. 5', G_5(x) ls 1 lw 2 notitle

set key top right
set yrange [0:0.05]
set xrange [0:45]
plot 'T2_v1.txt' ps 3 t 'SE, č. 1', F_1(x) ls 1 lw 2 notitle
set yrange [0:0.05]
set xrange [0:40]
plot 'T2_v2.txt' ps 3 t 'SE, č. 2', F_2(x) ls 1 lw 2 notitle
set yrange [0:0.05]
set xrange [0:40]
plot 'T2_v3.txt' ps 3 t 'SE, č. 3', F_3(x) ls 1 lw 2 notitle
set yrange [0:0.05]
set xrange [0:35]
plot 'T2_v4.txt' ps 3 t 'SE, č. 4', F_4(x) ls 1 lw 2 notitle
set yrange [0:0.035]
set xrange [0:25]
plot 'T2_v5.txt' ps 3 t 'SE, č. 5', F_5(x) ls 1 lw 2 notitle


unset multiplot
set terminal wxt
set output