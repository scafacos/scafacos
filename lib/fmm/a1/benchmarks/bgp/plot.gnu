set terminal postscript enhanced "Helvetica" 18 color

set size 1.0,1.0

########################################################

set key top left

set xlabel "dimension"
set ylabel "cycles"
set xrange [0:3000]
set yrange [0:30000]

set output "accumulate_unaligned_small.ps"

set title 'accumulate performance'

plot "190704.output" using 1:2 title 'aligned - generic' with lines, \
     "190704.output" using 1:3 title 'aligned - optimized' with lines, \
     "190704.output" using 1:4 title 'unaligned - generic' with lines, \
     "190704.output" using 1:5 title 'unaligned - optimized' with lines

########################################################

set key top left

set xlabel "dimension"
set ylabel "cycles"
set xrange [0:16384]
set yrange [0:125000]

set output "accumulate_unaligned_large.ps"

set title 'accumulate performance'

plot "190704.output" using 1:2 title 'aligned - generic' with lines, \
     "190704.output" using 1:3 title 'aligned - optimized' with lines, \
     "190704.output" using 1:4 title 'unaligned - generic' with lines, \
     "190704.output" using 1:5 title 'unaligned - optimized' with lines

########################################################
