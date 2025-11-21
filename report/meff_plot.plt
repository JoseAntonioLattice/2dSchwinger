reset 
set terminal epslatex color colortext standalone

set output "meff_plot.tex"

set size sq
set title "$V = 16^2$"
set xlabel "$t$"
set ylabel "$m_{\\tiny \\textrm{eff}}$"

set key Left top center reverse
p [-1:16]'effective_mass_m=0.dat' i 0 t "\\small $m_0 = 0$" lw 2, 'effective_mass_m=-0.05.dat' i 0 t "\\small $m_0 = -0.05$" lw 2

unset output
