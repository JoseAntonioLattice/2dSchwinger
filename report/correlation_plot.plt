reset 
set terminal epslatex color colortext 

set output "correlation_pion.tex"

L = 16
beta = 4

set size sq
set logscale y
set title sprintf("$V = %d^2, \\beta = %d$",L,beta)
set xlabel "$t$"
set ylabel "$C_{\\pi}(t)$"

set key Left top center reverse
f(x) = A*cosh(m*(x-L/2))
fit [4:11] f(x) '../data/pion_correlator_m=0.dat' u 1:2:3 via A, m
p [-1:17][0.004:0.2]'../data/pion_correlator_m=0.dat' w e t "\\small $m_0 = 0$" lw 2, [4:11] f(x) lw 2 t "\\small $A\\cosh(m_{\\pi}(t-L_t/2))$", keyentry t sprintf("\\small $m_{\\pi} = %1.3f(%d)$",m,1000*m_err)
unset output
