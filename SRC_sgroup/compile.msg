rm -f sgroup.o math_sg.o sto.o io.o lat.o pgrp.o pgrp_op.o  rotb.o sgrp_cub.o sgrp_hex.o sgrp_tet.o  sgrp_ort.o sgrp_trg.o sgrp_mon.o sgrp_trc.o 
icc       -c  sgroup.c
sgroup.c(1307): (col. 4) remark: LOOP WAS VECTORIZED.
sgroup.c(96): (col. 4) remark: LOOP WAS VECTORIZED.
sgroup.c(673): (col. 4) remark: LOOP WAS VECTORIZED.
sgroup.c(256): (col. 5) remark: LOOP WAS VECTORIZED.
sgroup.c(986): (col. 4) remark: LOOP WAS VECTORIZED.
icc       -c  math_sg.c
math_sg.c(276): (col. 4) remark: LOOP WAS VECTORIZED.
math_sg.c(140): (col. 4) remark: LOOP WAS VECTORIZED.
math_sg.c(115): (col. 3) remark: LOOP WAS VECTORIZED.
math_sg.c(107): (col. 3) remark: LOOP WAS VECTORIZED.
math_sg.c(56): (col. 4) remark: PARTIAL LOOP WAS VECTORIZED.
math_sg.c(58): (col. 4) remark: LOOP WAS VECTORIZED.
math_sg.c(44): (col. 4) remark: PARTIAL LOOP WAS VECTORIZED.
math_sg.c(46): (col. 7) remark: LOOP WAS VECTORIZED.
math_sg.c(17): (col. 4) remark: LOOP WAS VECTORIZED.
icc       -c  sto.c
icc       -c  io.c
io.c(1039): (col. 6) remark: LOOP WAS VECTORIZED.
icc       -c  lat.c
icc       -c  pgrp.c
icc       -c  pgrp_op.c
icc       -c  rotb.c
icc       -c  sgrp_cub.c
icc       -c  sgrp_hex.c
icc       -c  sgrp_tet.c
icc       -c  sgrp_ort.c
icc       -c  sgrp_trg.c
icc       -c  sgrp_mon.c
icc       -c  sgrp_trc.c
icc       -o sgroup sgroup.o math_sg.o sto.o io.o lat.o pgrp.o pgrp_op.o  rotb.o sgrp_cub.o sgrp_hex.o sgrp_tet.o  sgrp_ort.o sgrp_trg.o sgrp_mon.o sgrp_trc.o   -lm
