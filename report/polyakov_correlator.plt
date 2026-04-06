reset
#set terminal epslatex color colortext

#set output "correlation_polyakov.tex"

Lt = 8
Lx = 16
m0 = -0.1

#set title sprintf("$V = %d\\times%d, m_0 = %1.1f$",Lt,Lx,m0)
#set xlabel "$r$"
#set ylabel "$-\\log(C(r))/L_t$"

set title sprintf("V = %d×%d, m_0 = %1.1f",Lt,Lx,m0)
set xlabel "r"
set ylabel "-log(C(r))/L_t"


set key outside
#p for [j=1:9]'polyakov_correlator.dat' i j u 1:(-log($2)/8):(1/($2*8)*$3) w e t sprintf('$\\beta = %i$',j+1)
#unset output
p [-0.1:Lt][-0.01:1]for [j=1:9]'polyakov_correlator.dat' i j u 1:(-log($2)/Lt):(1/($2*Lt)*$3) w e t sprintf('{/Symbol b} = %i',j+1)
