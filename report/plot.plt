reset
set terminal epslatex color colortext standalone

set title "\\Large Schwinger model with two fermion flavors, $V = 8^2$"
set xlabel "$\\beta$"
set ylabel "plaquette value"
set key bottom Left title "$m_0$" enhanced left height 1.5

set output "plaquette.tex"

p [-0.1:4.1][0:1] 'm0_-0.1.dat' w e t "$-0.1$",\
'm0_0.0.dat' w e t "$0$", 'm0_10.dat' w e t "$10$",\
'm0_100.dat' w e t "$100$", 'm0_1000.dat' w e t "$1000$", 'quenched.dat' w e t "$\\infty$"
unset output
