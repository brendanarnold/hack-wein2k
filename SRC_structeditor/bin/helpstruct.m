 function helpstruct

printf("\nWien2k struct editor: LIST of FUNCTIONS:\n\n");
printf("a2adist			* calculates  distance between atoms\n");
printf("mina2adist		* calculates minimum distance between atoms\n");
printf("addatom			* adds an atom to the structure\n");
printf("addeqatom		* adds an atom and all equivalent\n");
printf("copyatom                * creates a copy of an atom\n");
printf("getaname		* converts atomic number into atomic symbol\n");
printf("getar0			* calculates r0 from atomic number\n");
printf("getazz			* converts atomic name into atomic number\n");
printf("loadstruct		* reads Wien2k structfile\n");
printf("movealla		* moves all atoms with vector vec\n");
printf("replaceatom		* replaces an atom with other atom\n");
printf("replaceeqatoms		* replaces an atom and all eqivalent with other atoms\n");
printf("rescale_c		* rescales z-positions for new c-lat.param. (vacuum at z=0.5)\n");
printf("rescale_c_2		* rescales z-positions for new c-lat.param. (vacuum at z=0.0)\n");
printf("rmatom			* removes an atom\n");
printf("rmeqatoms		* removes an atom and all equivalent\n");
printf("savestruct		* saves crystal structure\n");
printf("showequivalent		* outputs list of equivalent atoms\n");
printf("showstruct		* displays structure\n");
printf("smultatom		* creates symmetry equivalent positions\n");
printf("sshift			* symmetric shifts of equicalent atoms\n");
printf("makeconventional        * converts structure into the conventional form\n");
printf("makeprimitive           * converts structure to the primitive form\n");
printf("makesupercell           * creates supercell\n");
printf("makesurface             * creates surface for a given unitcells\n");
printf("\n");
printf("type \"help function_name\" for more help \n");

 end
