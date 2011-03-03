all:
	
	make -C SRC_ncmsymmetry 
	cp SRC_ncmsymmetry/ncmsymmetry bin/
	make -C SRC_readwrite
	cp SRC_readwrite/readwrite bin/
	make -C SRC_struct2mol
	cp SRC_struct2mol/struct2mol bin/
	make -C SRC_structgen 
	cp SRC_structgen/structgen bin/
clean:
	
	make clean -C SRC_ncmsymmetry 
	make clean -C SRC_readwrite
	make clean -C SRC_struct2mol
	make clean -C SRC_structgen 
