rm -f reformat.o 
icc -c  reformat.c
icc -o reformat reformat.o   -lm
