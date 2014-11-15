#set terminal wxt size 1400,700
set terminal png size 1800,900
# 1: yp, 2: heinz, 3: gal, 4: pim, 5: potato
# most frequent: col 3:          ((3,(4,(1,2))),5)
# 2nd most frequent: col 2: 	 ((4,(3,(1,2))),5)
# 3rd most frequent: col 4:     (((1,2),(3,4)),5)
# 4th most frequent: col 5:     ((2,(1,(3,4))),5)
# 5th most frequent: col 10:    ((1,(2,(3,4))),5)
# 6th most frequent: col 9:     (1,((2,(3,5)),4))
# 7th most frequent: col 7:     (1,((2,3),(4,5)))
# 8th most frequent: col 12:    (1,(((2,5),3),4))
# 9th most frequent: col 15:    (1,(((2,3),5),4))
# 10th most frequent: col 11:    (1,(((2,5),4),3))
# 11th most frequent: col 16:    (1,(((2,3),4),5))

 unset key
#set key bottom 
# set style data filledcurves x1
# plot [0:*] 'topo_ppd_vs_chunk1-120.out' using ($3+$4+$5+$6+$11), '' using ($3+$4+$5+$6), '' using ($3+$4+$5), '' using ($3+$4), '' using 4
set style data  histograms
set style histogram rowstacked
set style fill solid 1.0 border 1

# 1-840
set multiplot
set size 1,0.25
set origin 0.0, 0.75
plot [-0.5:839.5][0:6400] '1e4chunks_topo_pd_vs_chunks1-843' \
     using 3 t'((3,(4,(1,2))),5)', \
     '' using 2 t'((4,(3,(1,2))),5)', \
     '' using 4 t'(((1,2),(3,4)),5)', \
     '' using 5 t'((2,(1,(3,4))),5)', \
     '' using 10 t'((1,(2,(3,4))),5)', \
     '' using 9 t'(1,((2,(3,5)),4))', \
     '' using 7 t'(1,((2,3),(4,5)))', \
     '' using 12 t'(1,(((2,5),3),4))', \
     '' using 15 t'(1,(((2,3),5),4))', \
     '' using 11 t'(1,(((2,5),4),3))'
# pause -1

set origin 0.0, 0.5
plot [-0.5:419.5][0:6400] '2e4chunks_topo_pd_vs_chunks1-420' \
     using 3 t'((3,(4,(1,2))),5)', \
     '' using 2 t'((4,(3,(1,2))),5)', \
     '' using 4 t'(((1,2),(3,4)),5)', \
     '' using 5 t'((2,(1,(3,4))),5)', \
     '' using 10 t'((1,(2,(3,4))),5)', \
     '' using 9 t'(1,((2,(3,5)),4))', \
     '' using 7 t'(1,((2,3),(4,5)))', \
     '' using 12 t'(1,(((2,5),3),4))', \
     '' using 15 t'(1,(((2,3),5),4))', \
     '' using 11 t'(1,(((2,5),4),3))'
# pause -1

set origin 0.0, 0.25
plot [-0.5:209.5][0:6400] '4e4chunks_topo_pd_vs_chunks1-209' \
     using 3 t'((3,(4,(1,2))),5)', \
     '' using 2 t'((4,(3,(1,2))),5)', \
     '' using 4 t'(((1,2),(3,4)),5)', \
     '' using 5 t'((2,(1,(3,4))),5)', \
     '' using 10 t'((1,(2,(3,4))),5)', \
 '' using 9 t'(1,((2,(3,5)),4))', \
     '' using 7 t'(1,((2,3),(4,5)))', \
     '' using 12 t'(1,(((2,5),3),4))', \
     '' using 15 t'(1,(((2,3),5),4))', \
     '' using 11 t'(1,(((2,5),4),3))'

set origin 0.0, 0.0
plot [-0.5:104.5][0:6400] '8e4chunks_topo_pd_vs_chunks1-105' \
     using 3 t'((3,(4,(1,2))),5)', \
     '' using 2 t'((4,(3,(1,2))),5)', \
     '' using 4 t'(((1,2),(3,4)),5)', \
     '' using 5 t'((2,(1,(3,4))),5)', \
     '' using 10 t'((1,(2,(3,4))),5)', \
 '' using 9 t'(1,((2,(3,5)),4))', \
     '' using 7 t'(1,((2,3),(4,5)))', \
     '' using 12 t'(1,(((2,5),3),4))', \
     '' using 15 t'(1,(((2,3),5),4))', \
     '' using 11 t'(1,(((2,5),4),3))'

unset multiplot
pause -1


# 1-420
set multiplot
set size 1,0.25
set origin 0.0, 0.75
plot [-0.5:419.5][0:6400] '1e4chunks_topo_pd_vs_chunks1-843' \
     using 3 t'((3,(4,(1,2))),5)', \
     '' using 2 t'((4,(3,(1,2))),5)', \
     '' using 4 t'(((1,2),(3,4)),5)', \
     '' using 5 t'((2,(1,(3,4))),5)', \
     '' using 10 t'((1,(2,(3,4))),5)', \
     '' using 9 t'(1,((2,(3,5)),4))', \
     '' using 7 t'(1,((2,3),(4,5)))', \
     '' using 12 t'(1,(((2,5),3),4))', \
     '' using 15 t'(1,(((2,3),5),4))', \
     '' using 11 t'(1,(((2,5),4),3))'
# pause -1

set origin 0.0, 0.5
plot [-0.5:209.5][0:6400] '2e4chunks_topo_pd_vs_chunks1-420' \
     using 3 t'((3,(4,(1,2))),5)', \
     '' using 2 t'((4,(3,(1,2))),5)', \
     '' using 4 t'(((1,2),(3,4)),5)', \
     '' using 5 t'((2,(1,(3,4))),5)', \
     '' using 10 t'((1,(2,(3,4))),5)', \
     '' using 9 t'(1,((2,(3,5)),4))', \
     '' using 7 t'(1,((2,3),(4,5)))', \
     '' using 12 t'(1,(((2,5),3),4))', \
     '' using 15 t'(1,(((2,3),5),4))', \
     '' using 11 t'(1,(((2,5),4),3))'
# pause -1

set origin 0.0, 0.25
plot [-0.5:104.5][0:6400] '4e4chunks_topo_pd_vs_chunks1-209' \
     using 3 t'((3,(4,(1,2))),5)', \
     '' using 2 t'((4,(3,(1,2))),5)', \
     '' using 4 t'(((1,2),(3,4)),5)', \
     '' using 5 t'((2,(1,(3,4))),5)', \
     '' using 10 t'((1,(2,(3,4))),5)', \
 '' using 9 t'(1,((2,(3,5)),4))', \
     '' using 7 t'(1,((2,3),(4,5)))', \
     '' using 12 t'(1,(((2,5),3),4))', \
     '' using 15 t'(1,(((2,3),5),4))', \
     '' using 11 t'(1,(((2,5),4),3))'

set origin 0.0, 0.0
plot [-0.5:52][0:6400] '8e4chunks_topo_pd_vs_chunks1-105' \
     using 3 t'((3,(4,(1,2))),5)', \
     '' using 2 t'((4,(3,(1,2))),5)', \
     '' using 4 t'(((1,2),(3,4)),5)', \
     '' using 5 t'((2,(1,(3,4))),5)', \
     '' using 10 t'((1,(2,(3,4))),5)', \
 '' using 9 t'(1,((2,(3,5)),4))', \
     '' using 7 t'(1,((2,3),(4,5)))', \
     '' using 12 t'(1,(((2,5),3),4))', \
     '' using 15 t'(1,(((2,3),5),4))', \
     '' using 11 t'(1,(((2,5),4),3))'

unset multiplot
pause -1

# 420-840
set multiplot
set size 1,0.25
set origin 0.0, 0.75
plot [419.5:839.5][0:6400] '1e4chunks_topo_pd_vs_chunks1-843' \
     using 3 t'((3,(4,(1,2))),5)', \
     '' using 2 t'((4,(3,(1,2))),5)', \
     '' using 4 t'(((1,2),(3,4)),5)', \
     '' using 5 t'((2,(1,(3,4))),5)', \
     '' using 10 t'((1,(2,(3,4))),5)', \
     '' using 9 t'(1,((2,(3,5)),4))', \
     '' using 7 t'(1,((2,3),(4,5)))', \
     '' using 12 t'(1,(((2,5),3),4))', \
     '' using 15 t'(1,(((2,3),5),4))', \
     '' using 11 t'(1,(((2,5),4),3))'
# pause -1

set origin 0.0, 0.5
plot [209.5:419.5][0:6400] '2e4chunks_topo_pd_vs_chunks1-420' \
     using 3 t'((3,(4,(1,2))),5)', \
     '' using 2 t'((4,(3,(1,2))),5)', \
     '' using 4 t'(((1,2),(3,4)),5)', \
     '' using 5 t'((2,(1,(3,4))),5)', \
     '' using 10 t'((1,(2,(3,4))),5)', \
     '' using 9 t'(1,((2,(3,5)),4))', \
     '' using 7 t'(1,((2,3),(4,5)))', \
     '' using 12 t'(1,(((2,5),3),4))', \
     '' using 15 t'(1,(((2,3),5),4))', \
     '' using 11 t'(1,(((2,5),4),3))'
# pause -1

set origin 0.0, 0.25
plot [104.5:209.5][0:6400] '4e4chunks_topo_pd_vs_chunks1-209' \
     using 3 t'((3,(4,(1,2))),5)', \
     '' using 2 t'((4,(3,(1,2))),5)', \
     '' using 4 t'(((1,2),(3,4)),5)', \
     '' using 5 t'((2,(1,(3,4))),5)', \
     '' using 10 t'((1,(2,(3,4))),5)', \
 '' using 9 t'(1,((2,(3,5)),4))', \
     '' using 7 t'(1,((2,3),(4,5)))', \
     '' using 12 t'(1,(((2,5),3),4))', \
     '' using 15 t'(1,(((2,3),5),4))', \
     '' using 11 t'(1,(((2,5),4),3))'

set origin 0.0, 0.0
plot [52:104.5][0:6400] '8e4chunks_topo_pd_vs_chunks1-105' \
     using 3 t'((3,(4,(1,2))),5)', \
     '' using 2 t'((4,(3,(1,2))),5)', \
     '' using 4 t'(((1,2),(3,4)),5)', \
     '' using 5 t'((2,(1,(3,4))),5)', \
     '' using 10 t'((1,(2,(3,4))),5)', \
 '' using 9 t'(1,((2,(3,5)),4))', \
     '' using 7 t'(1,((2,3),(4,5)))', \
     '' using 12 t'(1,(((2,5),3),4))', \
     '' using 15 t'(1,(((2,3),5),4))', \
     '' using 11 t'(1,(((2,5),4),3))'

unset multiplot
pause -1





exit


set multiplot
set size 1,0.33
set origin 0.0, 0.67
plot [-0.5:21.5] '1e5chunks_topo_pd_vs_chunks1-84' \
     using 3 t'((3,(4,(1,2))),5)', \
     '' using 2 t'((4,(3,(1,2))),5)', \
     '' using 4 t'(((1,2),(3,4)),5)', \
     '' using 5 t'((2,(1,(3,4))),5)', \
     '' using 10 t'((1,(2,(3,4))),5)', \
'' using 9 t'(1,((2,(3,5)),4))', \
     '' using 7 t'(1,((2,3),(4,5)))', \
     '' using 12 t'(1,(((2,5),3),4))', \
     '' using 15 t'(1,(((2,3),5),4))', \
     '' using 11 t'(1,(((2,5),4),3))'
# pause -1

set origin 0.0, 0.33
plot [-0.5:52.5] '4e4chunks_topo_pd_vs_chunks1-209' \
     using 3 t'((3,(4,(1,2))),5)', \
     '' using 2 t'((4,(3,(1,2))),5)', \
     '' using 4 t'(((1,2),(3,4)),5)', \
     '' using 5 t'((2,(1,(3,4))),5)', \
     '' using 10 t'((1,(2,(3,4))),5)', \
 '' using 9 t'(1,((2,(3,5)),4))', \
     '' using 7 t'(1,((2,3),(4,5)))', \
     '' using 12 t'(1,(((2,5),3),4))', \
     '' using 15 t'(1,(((2,3),5),4))', \
     '' using 11 t'(1,(((2,5),4),3))'

set origin 0.0, 0.0
plot [-0.5:210.5] '1e4chunks_topo_pd_vs_chunks1-843' \
     using 3 t'((3,(4,(1,2))),5)', \
     '' using 2 t'((4,(3,(1,2))),5)', \
     '' using 4 t'(((1,2),(3,4)),5)', \
     '' using 5 t'((2,(1,(3,4))),5)', \
     '' using 10 t'((1,(2,(3,4))),5)', \
 '' using 9 t'(1,((2,(3,5)),4))', \
     '' using 7 t'(1,((2,3),(4,5)))', \
     '' using 12 t'(1,(((2,5),3),4))', \
     '' using 15 t'(1,(((2,3),5),4))', \
     '' using 11 t'(1,(((2,5),4),3))'
unset multiplot
pause -1



set multiplot
set size 1,0.33
set origin 0.0, 0.67
plot [20.5:42.5] '1e5chunks_topo_pd_vs_chunks1-84' \
     using 3 t'((3,(4,(1,2))),5)', \
     '' using 2 t'((4,(3,(1,2))),5)', \
     '' using 4 t'(((1,2),(3,4)),5)', \
     '' using 5 t'((2,(1,(3,4))),5)', \
     '' using 10 t'((1,(2,(3,4))),5)', \
'' using 9 t'(1,((2,(3,5)),4))', \
     '' using 7 t'(1,((2,3),(4,5)))', \
     '' using 12 t'(1,(((2,5),3),4))', \
     '' using 15 t'(1,(((2,3),5),4))', \
     '' using 11 t'(1,(((2,5),4),3))'

set origin 0.0, 0.33
plot [51.5:105.5] '4e4chunks_topo_pd_vs_chunks1-209' \
     using 3 t'((3,(4,(1,2))),5)', \
     '' using 2 t'((4,(3,(1,2))),5)', \
     '' using 4 t'(((1,2),(3,4)),5)', \
     '' using 5 t'((2,(1,(3,4))),5)', \
     '' using 10 t'((1,(2,(3,4))),5)', \
 '' using 9 t'(1,((2,(3,5)),4))', \
     '' using 7 t'(1,((2,3),(4,5)))', \
     '' using 12 t'(1,(((2,5),3),4))', \
     '' using 15 t'(1,(((2,3),5),4))', \
     '' using 11 t'(1,(((2,5),4),3))'

set origin 0.0, 0.0
plot [209.5:420.5] '1e4chunks_topo_pd_vs_chunks1-843' \
     using 3 t'((3,(4,(1,2))),5)', \
     '' using 2 t'((4,(3,(1,2))),5)', \
     '' using 4 t'(((1,2),(3,4)),5)', \
     '' using 5 t'((2,(1,(3,4))),5)', \
     '' using 10 t'((1,(2,(3,4))),5)', \
'' using 9 t'(1,((2,(3,5)),4))', \
     '' using 7 t'(1,((2,3),(4,5)))', \
     '' using 12 t'(1,(((2,5),3),4))', \
     '' using 15 t'(1,(((2,3),5),4))', \
     '' using 11 t'(1,(((2,5),4),3))'
unset multiplot
pause -1

set multiplot
set size 1,0.33
set origin 0.0, 0.67
plot [41.5:63.5] '1e5chunks_topo_pd_vs_chunks1-84' \
     using 3 t'((3,(4,(1,2))),5)', \
     '' using 2 t'((4,(3,(1,2))),5)', \
     '' using 4 t'(((1,2),(3,4)),5)', \
     '' using 5 t'((2,(1,(3,4))),5)', \
     '' using 10 t'((1,(2,(3,4))),5)', \
  '' using 9 t'(1,((2,(3,5)),4))', \
     '' using 7 t'(1,((2,3),(4,5)))', \
     '' using 12 t'(1,(((2,5),3),4))', \
     '' using 15 t'(1,(((2,3),5),4))', \
     '' using 11 t'(1,(((2,5),4),3))'
# pause -1

set origin 0.0, 0.33
plot [104.5:157.5] '4e4chunks_topo_pd_vs_chunks1-209' \
     using 3 t'((3,(4,(1,2))),5)', \
     '' using 2 t'((4,(3,(1,2))),5)', \
     '' using 4 t'(((1,2),(3,4)),5)', \
     '' using 5 t'((2,(1,(3,4))),5)', \
     '' using 10 t'((1,(2,(3,4))),5)', \
 '' using 9 t'(1,((2,(3,5)),4))', \
     '' using 7 t'(1,((2,3),(4,5)))', \
     '' using 12 t'(1,(((2,5),3),4))', \
     '' using 15 t'(1,(((2,3),5),4))', \
     '' using 11 t'(1,(((2,5),4),3))'

set origin 0.0, 0.0
plot [419.5:630.5] '1e4chunks_topo_pd_vs_chunks1-843' \
     using 3 t'((3,(4,(1,2))),5)', \
     '' using 2 t'((4,(3,(1,2))),5)', \
     '' using 4 t'(((1,2),(3,4)),5)', \
     '' using 5 t'((2,(1,(3,4))),5)', \
     '' using 10 t'((1,(2,(3,4))),5)', \
 '' using 9 t'(1,((2,(3,5)),4))', \
     '' using 7 t'(1,((2,3),(4,5)))', \
     '' using 12 t'(1,(((2,5),3),4))', \
     '' using 15 t'(1,(((2,3),5),4))', \
     '' using 11 t'(1,(((2,5),4),3))'

unset multiplot
pause -1

set multiplot
set size 1,0.33
set origin 0.0, 0.67
plot [62.5:84.5] '1e5chunks_topo_pd_vs_chunks1-84' \
     using 3 t'((3,(4,(1,2))),5)', \
     '' using 2 t'((4,(3,(1,2))),5)', \
     '' using 4 t'(((1,2),(3,4)),5)', \
     '' using 5 t'((2,(1,(3,4))),5)', \
     '' using 10 t'((1,(2,(3,4))),5)', \
 '' using 9 t'(1,((2,(3,5)),4))', \
     '' using 7 t'(1,((2,3),(4,5)))', \
     '' using 12 t'(1,(((2,5),3),4))', \
     '' using 15 t'(1,(((2,3),5),4))', \
     '' using 11 t'(1,(((2,5),4),3))'
# pause -1

set origin 0.0, 0.33
plot [156.5:210.5] '4e4chunks_topo_pd_vs_chunks1-209' \
     using 3 t'((3,(4,(1,2))),5)', \
     '' using 2 t'((4,(3,(1,2))),5)', \
     '' using 4 t'(((1,2),(3,4)),5)', \
     '' using 5 t'((2,(1,(3,4))),5)', \
     '' using 10 t'((1,(2,(3,4))),5)', \
 '' using 9 t'(1,((2,(3,5)),4))', \
     '' using 7 t'(1,((2,3),(4,5)))', \
     '' using 12 t'(1,(((2,5),3),4))', \
     '' using 15 t'(1,(((2,3),5),4))', \
     '' using 11 t'(1,(((2,5),4),3))'


set origin 0.0, 0.0
plot [629.5:843.5] '1e4chunks_topo_pd_vs_chunks1-843' \
     using 3 t'((3,(4,(1,2))),5)', \
     '' using 2 t'((4,(3,(1,2))),5)', \
     '' using 4 t'(((1,2),(3,4)),5)', \
     '' using 5 t'((2,(1,(3,4))),5)', \
     '' using 10 t'((1,(2,(3,4))),5)', \
 '' using 9 t'(1,((2,(3,5)),4))', \
     '' using 7 t'(1,((2,3),(4,5)))', \
     '' using 12 t'(1,(((2,5),3),4))', \
     '' using 15 t'(1,(((2,3),5),4))', \
     '' using 11 t'(1,(((2,5),4),3))'

unset multiplot

pause -1

exit

