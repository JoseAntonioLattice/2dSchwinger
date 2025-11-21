reset 
set terminal epslatex color colortext standalone

set output "correlation_pion.tex"

set size sq
set logscale y
set title "$V = 16^2$"
set xlabel "$t$"
set ylabel "pion correlation"

set key Left top center reverse
p [-1:17][0.004:0.2]'pion_correlator_m=0.dat' w e t "\\small $m_0 = 0$" lw 2, 'pion_correlator.dat' w e t '\small $m_0 = -0.05$' lw 2
unset output
