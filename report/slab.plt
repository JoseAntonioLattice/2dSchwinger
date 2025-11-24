reset

set terminal epslatex color colortext standalone

beta = 3
L = 16
V = L**2
Q = 0
file = "../data/topological_charge_1.dat"
set output sprintf("top_char_slab_%i.tex",beta)

f(x) = chi*V*x*(1-x) - Q**2*x**2

set title sprintf("$Q = 0$, $\\beta = %d$",beta)
set xlabel "$x$"
set ylabel "$\\langle q^2 \\rangle$"

set key bottom

fit [0:.9999] f(x) file i 0 u ($1/L):2:3 via chi
p [0:1.01] file i 0 u ($1/L):2:3 w e notitle ls 4 lc 6 lw 2, keyentry t sprintf("$\\chi_t = %1.4f$",chi), f(x) t "$\\chi_t V x (1-x)$" lc 4 lw 2
