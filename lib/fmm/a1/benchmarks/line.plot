#!/usr/bin/gnuplot
set terminal postscript enhanced "Helvetica" 20
set pointsize 2
set data style lp
set key left top
set logscale x 2
set term postscript
set grid xtics ytics mxtics mytics
set style line 1 lt 1 lw 2
set style line 2 lt 2 lw 2
set style line 3 lt 3 lw 2
set style line 4 lt 4 lw 2
set style line 5 lt 5 lw 2
set style line 6 lt 6 lw 2

set yrange [0:]
set xrange [8:2097152]
set xtics (8,64,512,"4K" 4096,"32K" 32768,"256K" 262144,"2M" 2097152)
set xlabel 'Message size (bytes)'
set ylabel 'ARMCI\_Acc Bandwidth (MBPS)'
set output './ARMCI_Progress_mode.eps'
plot './armci/Performance/ARMCI_Accumulate_bw.data' using 1:4 index 0 title 'none' with lp lt 1 lw 5 pt 1, \
     './armci/Performance/ARMCI_Accumulate_bw.data' using 1:2 index 0 title 'interrupts' with lp lt 3 lw 5 pt 2, \
     './armci/Performance/ARMCI_Accumulate_bw.data' using 1:3 index 0 title 'helper thread' with lp lt 4 lw 5 pt 3

set yrange [0:]
set xrange [8:16384]
set xlabel 'Message size (bytes)'
set ylabel 'Get Latency (usec)'
set output './ARMCI_Get_latency.eps'
set xtics (8,16,32,64,128,256,512,1024,2048,4096,8192,16384)
plot './armci/Performance/ARMCI_Get_latency_comparision.data' using 1:2 index 0 title 'DCMF' with lp lt 1 lw 5, \
     './armci/Performance/ARMCI_Get_latency_comparision.data' using 1:5 index 0 title 'ARMCI (active target)' with lp lt 3 lw 5, \
     './armci/Performance/ARMCI_Get_latency_comparision.data' using 1:3 index 0 title 'ARMCI (interrupts)' with lp lt 4 lw 5, \
     './armci/Performance/ARMCI_Get_latency_comparision.data' using 1:4 index 0 title 'ARMCI (helper thread)' with lp lt 2 lw 5

set yrange [:]
set xrange [4096:262144]
set xtics ("512x1" 4096, "512x2" 8192, "512x4" 16384, "512x8" 32768, "512x16" 65536, "512x32" 131072, "512x64" 262144)
set xlabel 'Sub-matrix size (doubles)'
set ylabel 'Total Time for 8NbPutS + (DGEMM) + Wait in Cycles'
set output './A1_PutS_handoff.eps'

plot './a1/Performance/A1_PutS_handoff_512.data' using 1:5 index 0 title 'No-Handoff-Base' with lp lt 1 lw 5 pt 1, \
     './a1/Performance/A1_PutS_handoff_512.data' using 1:6 index 0 title 'No-Handoff+DGEMM-Overlap' with lp lt 3 lw 5 pt 2, \
     './a1/Performance/A1_PutS_handoff_512.data' using 1:3 index 0 title 'Handoff-Base' with lp lt 4 lw 5 pt 3, \
     './a1/Performance/A1_PutS_handoff_512.data' using 1:4 index 0 title 'Handoff+DGEMM-Overlap' with lp lt 2 lw 5 pt 4

set yrange [0:]
set xrange [8:32768]
set xtics (8,16,32,64,128,256,512,"1K" 1024,"2K" 2048,"4K" 4096,"8K" 8192,"16K" 16384,"32K" 32768)
set xlabel 'Message size (bytes)'
set ylabel 'Get Latency (usec)'
set output './Get_latency.eps'
plot './armci-vs-a1/Get_latency.data' using 1:2 index 0 title 'DCMF' with lp lt 1 lw 5 pt 1, \
     './armci-vs-a1/Get_latency.data' using 1:3 index 0 title 'A1' with lp lt 3 lw 5 pt 3, \
     './armci-vs-a1/Get_latency.data' using 1:4 index 0 title 'ARMCI' with lp lt 5 lw 5 pt 2

set yrange [0:]
set xrange [256:65536]
set xtics (256,512,"1K" 1024,"2K" 2048,"4K" 4096,"8K" 8192,"16K" 16384,"32K" 32768, "64K" 65536)
set xlabel 'Message size (bytes)'
set ylabel 'Get Bandwidth (MBPS)'
set output './Get_bw.eps'
plot './armci-vs-a1/Get_bw.data' using 1:2 index 0 title 'A1' with lp lt 3 lw 5 pt 3, \
     './armci-vs-a1/Get_bw.data' using 1:3 index 0 title 'ARMCI' with lp lt 5 lw 5 pt 2

set yrange [0:]
set xrange [4096:524288]
set xtics ("512x1" 4096, "512x2" 8192, "512x4" 16384, "512x8" 32768, "512x16" 65536, "512x32" 131072, "512x64" 262144, "512x128" 524288) 
set xlabel 'Sub-matrix size (doubles)'
set ylabel 'GetS Latency (usec)'
set output './GetS_latency_512.eps'
plot './armci-vs-a1/GetS_latency_512.data' using 1:3 index 0 title 'A1' with lp lt 3 lw 5 pt 3, \
     './armci-vs-a1/GetS_latency_512.data' using 1:4 index 0 title 'ARMCI' with lp lt 5 lw 5 pt 2

set yrange [0:]
set xrange [128:524288]
set xtics ("4x4" 128, "8x8" 512, "16x16" 2048, "32x32" 8192, "64x64" 32768, "128x128" 131072, "256x256" 524288)
set xlabel 'Sub-matrix size (doubles)'
set ylabel 'GetS Bandwidth (MBPS)'
set output './GetS_bw.eps'
plot './armci-vs-a1/GetS_bw.data' using 1:3 index 0 title 'A1' with lp lt 3 lw 5 pt 3, \
     './armci-vs-a1/GetS_bw.data' using 1:4 index 0 title 'ARMCI' with lp lt 5 lw 5 pt 2

set yrange [0:]
set xrange [4096:1048576]
set xtics ("512x1" 4096, "512x4" 16384, "512x16" 65536, "512x64" 262144, "512x256" 1048576)
set xlabel 'Sub-matrix size (doubles)'
set ylabel 'GetS Bandwidth (MBPS)'
set output './GetS_bw_512.eps'
plot './armci-vs-a1/GetS_bw_512.data' using 1:3 index 0 title 'A1' with lp lt 3 lw 5 pt 3, \
     './armci-vs-a1/GetS_bw_512.data' using 1:4 index 0 title 'ARMCI' with lp lt 5 lw 5 pt 2

set yrange [0:]
set xrange [8:32768]
set xtics (8,16,32,64,128,256,512,"1K" 1024,"2K" 2048,"4K" 4096,"8K" 8192,"16K" 16384,"32K" 32768)
set xlabel 'Message size (bytes)'
set ylabel 'Put Latency - local completion (usec)'
set output './Put_latency_local.eps'
plot './armci-vs-a1/Put_latency_local.data' using 1:2 index 0 title 'DCMF' with lp lt 1 lw 5 pt 1, \
     './armci-vs-a1/Put_latency_local.data' using 1:3 index 0 title 'A1' with lp lt 3 lw 5 pt 3, \
     './armci-vs-a1/Put_latency_local.data' using 1:4 index 0 title 'ARMCI' with lp lt 5 lw 5 pt 2

set yrange [0:]
set xrange [8:32768]
set xtics (8,16,32,64,128,256,512,"1K" 1024,"2K" 2048,"4K" 4096,"8K" 8192,"16K" 16384,"32K" 32768)
set xlabel 'Message size (bytes)'
set ylabel 'Put Latency - remote completion (usec)'
set output './Put_latency_remote.eps'
plot './armci-vs-a1/Put_latency_remote.data' using 1:2 index 0 title 'DCMF' with lp lt 1 lw 5 pt 1, \
     './armci-vs-a1/Put_latency_remote.data' using 1:3 index 0 title 'A1' with lp lt 3 lw 5 pt 3, \
     './armci-vs-a1/Put_latency_remote.data' using 1:4 index 0 title 'ARMCI' with lp lt 5 lw 5 pt 2

set yrange [0:]
set xrange [256:65536]
set xtics (256,512,"1K" 1024,"2K" 2048,"4K" 4096,"8K" 8192,"16K" 16384,"32K" 32768, "64K" 65536)
set xlabel 'Message size (bytes)'
set ylabel 'Put Bandwidth (MBPS)'
set output './Put_bw.eps'
plot './armci-vs-a1/Put_bw.data' using 1:2 index 0 title 'A1' with lp lt 3 lw 5 pt 3, \
     './armci-vs-a1/Put_bw.data' using 1:3 index 0 title 'ARMCI' with lp lt 5 lw 5 pt 2


set yrange [0:]
set xrange [4096:524288]
set xtics ("512x1" 4096, "512x2" 8192, "512x4" 16384, "512x8" 32768, "512x16" 65536, "512x32" 131072, "512x64" 262144, "512x128" 524288)
set xlabel 'Sub-Matrix size (doubles)'
set ylabel 'AccS Local Completion Latency (usec)'
set output './AccS_latency_local.eps'
plot './armci-vs-a1/PutAccS_latency.data' using 1:3 index 0 title 'A1' with lp lt 3 lw 5 pt 3, \
     './armci-vs-a1/PutAccS_latency.data' using 1:5 index 0 title 'ARMCI' with lp lt 5 lw 5 pt 2

set yrange [0:]
set xrange [4096:524288]
set xtics ("512x1" 4096, "512x2" 8192, "512x4" 16384, "512x8" 32768, "512x16" 65536, "512x32" 131072, "512x64" 262144, "512x128" 524288)
set xlabel 'Sub-matrix size (doubles)'
set ylabel 'AccS Remote Completion Latency (usec)'
set output './AccS_latency_remote.eps'
plot './armci-vs-a1/PutAccS_latency.data' using 1:4 index 0 title 'A1' with lp lt 3 lw 5 pt 3, \
     './armci-vs-a1/PutAccS_latency.data' using 1:6 index 0 title 'ARMCI' with lp lt 5 lw 5 pt 2
