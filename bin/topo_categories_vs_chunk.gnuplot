set terminal wxt size 1600,800
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

# instead of treating each topology as a separate category
# lump them together in a meaningful way:

# 1) categorize them according to the taxon (or taxa) in the other child
# of heinz's parent in the tree (rooted with potato as outgroup)
# 1,2,3 (cols 2,3,4)   this is 1 (yellow pear)  "2-1"
# 5,12,13 this is 4 (pimpinellifolium) "2-4"
# 6,14,15 this is 3 (galapagense) "2-3"
# 4,7,8,9,10,11 this is a tree containing multiple taxa (1,3,4 or two of those 3) "2-multi" 

# 2) categorize them according to which species gal or pim (3 or 4) has the more recent
# last common ancestor with yellow pear (1).
# galapagense (3):        1,6,7,10,12;
# pimpinellifolium (4):   2,5,8,11,14;
# gal and pim have the same last common ancestor with yellow pear: 3,4,9,13,15.

# unset key
set key bottom 
# set style data filledcurves x1
# plot [0:*] 'topo_ppd_vs_chunk1-120.out' using ($3+$4+$5+$6+$11), '' using ($3+$4+$5+$6), '' using ($3+$4+$5), '' using ($3+$4), '' using 4
set style data  histograms
set style histogram rowstacked
set style fill solid 1.0 border 1


# find most recent common ancestor
# of 1 and 3, 1 and 4;
# which is most recent?
# possible answers are 3, 4, equal.
# 1-840
set multiplot
set size 1,0.25
set origin 0.0, 0.75
plot [-0.5:839.5][0:6400] '1e4chunks_topo_pd_vs_chunks1-843' \
     using ($2+$7+$8+$11+$13) t'(yp,gal)', \
     '' using ($3+$6+$9+$12+$15) t'(yp,pim)', \
     '' using ($4+$5+$10+$14+$16) t'(yp,(gal,pim))'

set origin 0.0, 0.5
plot [-0.5:419.5][0:6400] '2e4chunks_topo_pd_vs_chunks1-420' \
      using ($2+$7+$8+$11+$13) t'(yp,gal)', \
     '' using ($3+$6+$9+$12+$15) t'(yp,pim)', \
     '' using ($4+$5+$10+$14+$16) t'(yp,(gal,pim))'

# pause -1

set origin 0.0, 0.25
plot [-0.5:209.5][0:6400] '4e4chunks_topo_pd_vs_chunks1-209' \
   using ($2+$7+$8+$11+$13) t'(yp,gal)', \
     '' using ($3+$6+$9+$12+$15) t'(yp,pim)', \
     '' using ($4+$5+$10+$14+$16) t'(yp,(gal,pim))'
  
set origin 0.0, 0.0
plot [-0.5:104.5][0:6400] '8e4chunks_topo_pd_vs_chunks1-105' \
  using ($2+$7+$8+$11+$13) t'(yp,gal)', \
     '' using ($3+$6+$9+$12+$15) t'(yp,pim)', \
     '' using ($4+$5+$10+$14+$16) t'(yp,(gal,pim))'
  
unset multiplot
pause -1

# which
# 1-840
set multiplot
set size 1,0.25
set origin 0.0, 0.75
plot [-0.5:839.5][0:6400] '1e4chunks_topo_pd_vs_chunks1-843' \
     using ($3+$4+$2) t'(2,1)', \
     '' using ($6+$13+$14) t'(2,4)', \
     '' using ($7+$15+$16) t'(2,3)', \
     '' using ($5+$8+$9+$10+$11+$12) t'2,()'

set origin 0.0, 0.5
plot [-0.5:419.5][0:6400] '2e4chunks_topo_pd_vs_chunks1-420' \
     using ($3+$4+$2) t'(2,1)', \
     '' using ($6+$13+$14) t'(2,4)', \
     '' using ($7+$15+$16) t'(2,3)', \
     '' using ($5+$8+$9+$10+$11+$12) t'2,()'
# pause -1

set origin 0.0, 0.25
plot [-0.5:209.5][0:6400] '4e4chunks_topo_pd_vs_chunks1-209' \
  using ($3+$4+$2) t'(2,1)', \
     '' using ($6+$13+$14) t'(2,4)', \
     '' using ($7+$15+$16) t'(2,3)', \
     '' using ($5+$8+$9+$10+$11+$12) t'2,()'
  
set origin 0.0, 0.0
plot [-0.5:104.5][0:6400] '8e4chunks_topo_pd_vs_chunks1-105' \
using ($3+$4+$2) t'(2,1)', \
     '' using ($6+$13+$14) t'(2,4)', \
     '' using ($7+$15+$16) t'(2,3)', \
     '' using ($5+$8+$9+$10+$11+$12) t'2,()'
   
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

