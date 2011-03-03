rm -f module.o makeP1.o convert.o readtext.o 
clean
ifort -c -FR -mp -w module.f
ifort -c -FR -mp -w makeP1.f
ifort -c -FR -mp -w convert.f
ifort -c -FR -mp -w readtext.f
ifort -o ./makeP1  -L/area51/opt/intel/mkl721/lib/32 -static-libcxa      module.o makeP1.o convert.o readtext.o
