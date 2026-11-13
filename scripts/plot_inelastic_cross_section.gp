set title "Elastic electron scattering (PWBA)"
set xlabel "Scattering angle theta (deg)"
set ylabel "d sigma / d Omega (mb/sr)"

set logscale y
set grid
set key top right

plot "cross_section.dat" using 1:4 \
    with lines linewidth 2 \
    title "Uniform spherical charge distribution"
