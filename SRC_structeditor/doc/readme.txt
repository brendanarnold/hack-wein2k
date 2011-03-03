Install
-------

*Be sure octave and DataExplorer (DX) are installed (either from your 
Linux distribution or http://www.octave.org and http://www.opendx.org ).
*Compile fortran sources (SRC_*) and copy binaries to bin.
 To do this adjust options in make.inc file and run install script.
*Make sure that octave will find the bin directory, e.g. update
 octave enviroment variables:  OCTAVE_EXEC_PATH, OCTAVE_PATH.

setenv STRUCTEDIT_PATH $WIENROOT/SRC_structeditor/bin
set path = ($WIENROOT $STRUCTEDIT_PATH $path .)
setenv OCTAVE_EXEC_PATH ${PATH}::
setenv OCTAVE_PATH ${STRUCTEDIT_PATH}::
 

Usage
-----

The package is a collection of octave programs plus few Fortran
programs and shell scripts that are intended to provide WIEN2k struct
file editing and visualization possibility. It is NOT intended for the
novel (unexperienced) WIEN2k user since it requires both, a good idea 
of what you actually want to model and how WIEN2k stores it's structure, 
but the octave environment ensures full and migthy functionality and 
flexibility in terms of three main tasks any structure manipulation needs:

* loading the structure file and storing it as a convenient 
  data structure
* manipulation of this structure (e.g. adding, removing replacing 
  atoms, creating super-cells, surfaces, etc.) 
* visualization and saving structure in a form of Wien2k struct file.

Below these three tasks are in details described.

1. Loading and structure of structure.

After starting octave (use simply    octave   at the command line)
the Wien2k struct file is loaded with loadstruct function into 
an octave variable. 

  s=loadstruct("GaN.struct")
  s1=loadstruct("rh.struct");

The ";" at the end of the line restricts the output of a command

Here Gan.struct file is loaded into variable s (structure variable),
and "rh.struct" into variable s1. The structure variable (s or s1)
contains most important data describing unit cell:

 a        -  unit cell parameters , 
 alpha    -  lattice angles,
 brlat    -  Bravais lattice (row is a vector convention),
 lat2car  -  lattice to Cartesian transformation matrix (in  
             some cases it is different then Bravais matrix - 
             row is a vector convention),
 lattic   -  lattice type, 
 pos      -  positions of atoms, (all atoms in the structure, the
             concept of "equivalent" atoms is not used here) 
 aname    -  name of atoms, 
 zz       -  atomic numbers, 
 rmt      -  lapw sphere radii,
 r0       -  first point in radial mesh,
 jrj      -   number of points in radial mesh.

Example of GaN.struct stored as s.

s =
{
  a =

    6.02822700000000  6.02822700000000  9.80573700000000

  alpha =

     90   90  120

  aname =

Ga
Ga
N
N

  brlat =

     5.220597722000000  -3.014113500000000   0.000000000000000
     0.000000000000000   6.028227000000000   0.000000000000000
     0.000000000000000   0.000000000000000   9.805737000000001

  jrj =

    781
    781
    781
    781

  lat2car =

     5.220597722000000  -3.014113500000000   0.000000000000000
     0.000000000000000   6.028227000000000   0.000000000000000
     0.000000000000000   0.000000000000000   9.805737000000001

  lattic = H
  nat = 4
  pos =

    0.333333330000000  0.666666660000000  0.500000000000000
    0.666666660000000  0.333333330000000  0.000000000000000
    0.333333330000000  0.666666660000000  0.857000000000000
    0.666666660000000  0.333333330000000  0.357000000000000

  r0 =

     5.00000000000000e-05
     5.00000000000000e-05
     1.00000000000000e-04
     1.00000000000000e-04

  rmt =

    1.95000000000000
    1.95000000000000
    1.50000000000000
    1.50000000000000

  zz =

    31
    31
     7
     7
}

Because whole structure is stored as a single variable, it can be 
copied to another variable with simple assignment operation:

  s1=s

here s1 contains exactly the save data as s. This may be useful for
making temporal backup of introduced changes. Each part of the data
contained in the structure variables can be directly accessed using
usual structure oriented syntax, e.g.:

s.pos(3,1)=s.pos(3,1) + 0.1
aname(3)="bleble"

first coordinate of the third atom is increased by 0.1, and its name 
is changed to "bleble". Octave environment makes very simple any 
manipulation concerning vector and matrixes:

vec(:)=s.pos(3,:)
vec(:)=lat2car*vec

here position of the third atom is stored into vector vec, and this
vector is converted to Cartesian coordinates.

In order to simplify most common actions a number of functions has
been implemented. Here they are listed and shortly described, more
extended description can be found at the end of this document.

a2adist			* calculates  distance between atoms        
mina2adist		* calculates minimum distance between atoms
addatom			* adds an atom to the structure
addeqatom		* adds an atom and all equivalent
copyatom                * creates a copy of an atom
getaname		* converts atomic number into atomic symbol
getar0			* calculates r0 from atomic number
getazz			* converts atomic name into atomic number
loadstruct		* reads Wien2k structfile
movealla		* moves all atoms with vector vec
replaceatom		* replaces an atom with other atom
replaceeqatoms		* replaces an atom and all equivalent with other atoms
rescale_c               * rescales z-positions for new c-lat.param. (vacuum at z=0.5)
rescale_c_2             * rescales z-positions for new c-lat.param. (vacuum at z=0.0)
rmatom			* removes an atom
rmeqatoms		* removes an atom and all equivalent
savestruct		* saves crystal structure
showequivalent		* outputs list of equivalent atoms
showstruct		* displays structure (using DX)
smultatom		* creates symmetry equivalent positions
sshift			* symmetric shifts of equivalent atoms
makeconventional        * converts structure into the conventional form
makeprimitive           * converts structure to the primitive form
makesupercell           * creates supercell
makesurface             * creates surface for a given unitcell

Most of these functions perform some action with an input structure and
return a modified structure. They are usually called with a syntax:

       sout=functionname(sin,....),

where sout is an output structure, sin is an input structure, usually
followed by other arguments. Functions like: a2adist, mina2adist, 
getaname, getar0, getazz return only scalar or string variable.
Function showstruct does not return anything. 

Structures can be displayed with function showstruct, by simply typing:

      showstruct(sin),

where sin is an input structure. showstruct executes DataExplorer, thus 
it must be installed and accessible by octave. 

Structures can be saved into Wien2k file with function savestruct:

      savestruct(sin,"GaN.struct").

All functions are documented, using octave style (first comented out
lines are displayed by help engine). Type:

helpstruct

to get list of functions, and type:

help functionname

to get more detailed info about functionname.

A few examples that make use of above listed functions.

a) 

# load GaN wurtzite 
#
s=loadstruct ("GaN.struct");

# add He (for fun)
#
s1=addatom (s,"He",[0.1 0.1 0.1],0);

# fill also symmetry equivalent positions
#
s2=smultatom (s1,5);

# shift one He atom with [0.2 0.2 0.2] and other equivalent
# with equivalent directions
#
s3=sshift (s2,16,[0.2 0.2 0.2]);

# replace one He atom and all equivalent with H
#
s4=replaceeqatoms (s3,5,"H");

# display it
#
showstruct (s4);

# rm one atom and all equivalent
#
s5=rmeqatoms (s4,5); 

# change unitcell
a=[1 0 0; 1 1 0; 0 0 1]
s6=makesupercell (s4,a);

# save it as test.struct
#
savestruct (s6,"test.struct");

b) 

#load Rh fcc unitcell, 
#
s=loadstruct ("rh.struct");

#create 111 surface with 26 bohr "thickness" and 20 bohr vacuum, 
#makesurface works only for noncentred lattice,
#thus before that convert it to primitive
#
s1=makeprimitive (s);
n=[1 1 1];
s2=makesurface (s1,n,0,26,20);

# shift all atoms (into center of unit cell) for better view
#
shift=(s2.a(3)-a2adist (s2,7,1))/2
vec=[0 ,0 ,shift/s2.a(3)]
s3=movealla (s2,vec);


# add BN molecule on top (with 3 bohr Rh-N distance), 
#
pos=[0 0 s3.pos(1,3)-3/s3.a(3)];
s4=addeqatom (s3,"N",pos,0);
pos=[1/3 1/3 s3.pos(1,3)-3/s3.a(3)];
s5=addeqatom (s4,"B",pos,0);

# and create 3x3x1 supercel out of it.
#
a=[3 0 0; 0 3 0; 0 0 1];
s6=makesupercell (s5,a);

showstruct (s6);
c)

#lets remove one B (atom 95) from the top of surface created in the previous example
#and put three H atoms to saturate broken BN bonds
#
pos1=pos+(s6.pos(77,:)-s6.pos(95,:))*0.5
pos2=pos+(s6.pos(80,:)-s6.pos(95,:))*0.5
pos3=pos+(s6.pos(78,:)-s6.pos(95,:))*0.5
s7=addatom (s6,"H",pos1,0);
s7=addatom (s7,"H",pos2,0);
s7=addatom (s7,"H",pos3,0);
s7=rmatom (s7,95);

# or better to save inversion symmetry instead of the last 4 lines:
s7=addeqatom (s6,"H",pos1,0,"P-1");
s7=addeqatom (s7,"H",pos2,0,"P-1");
s7=addeqatom (s7,"H",pos3,0,"P-1");
s7=rmeqatoms (s7,95,"P-1");

 
Descriptions of implemented functions
-------------------------------------

   usage:   dis=a2adist(s,ind1,ind2)

         calculates  distance between atoms ind1 and ind2,
         periodic images are not taken into account

         dis       result
         s         input structure
         ind1      index of first atom
         ind2      index of second atom

   example:

       dis=a2adist(s,1,2)

******************************************************************

   usage:   dis=mina2adist(s,ind1,ind2)

         calculates minimum distance between atoms ind1 and ind2 
         taking into account periodic images         
 
         dis       result
         s         input structure
         ind1      index of first atom
         ind2      index of second atom

   example:    

       dis=mina2adist(s,1,2)  

*******************************************************************

  usage:  sr=addatom(s,atom,pos,ind) 

        adds an atom to the structure s and returns all as structure sr

        sr      output structure
        s       matlab structure used to store structural information
        atom    identifies atom, can be string variable containing
                atomic symbol or scalar with atomic number
        pos     position in unit cell
        ind     a index of new atom (position in the list),
                if 0 then atom is added to the end of the list

  example:
  
         sr=addatom(s,"N",[0,0,0.1],3)        

*******************************************************************

  usage:  sr=addeqatom(s,atom,pos,ind)

        adds an atom and all equivalent to the structure s and returns all as structure sr

        sr      output structure
        s       matlab structure used to store structural information
        atom    identifies atom, can be string variable containing
                atomic symbol or scalar with atomic number
        pos     position in unit cell
        ind     a index of new atom (position in the list),
                if 0 then atom is added to the end of the list
        spgn    spacegroup name, if not specified sym. operations
                are generated         
 
  example:

         sr=addeqatom(s,"N",[0,0,0.1],3,"P63mc")

*******************************************************************

   usage:  sr=copyatom(s,ind1,ind2,pos)

        creates a copy of atom ind1 and puts it in the structure with
        index ind2

        s     input structure
        ind1  index of original atom
        ind2  index of new atom (if = 0, then adds it at the end)
        pos   position of new atom

   example: sr=copyatom(s,1,0,[0.1 0.1 0.1})

*******************************************************************

   usage:  elem=getaname(elmzz)

          converts atomic number elemzz into atomic symbol elem

   example:  

          elem=getaname(5)

*******************************************************************

     usage:  r0=getar0(zz)  
    
          calculates r0 from atomic number zz

     example:

           r0=getar0(5)

*******************************************************************

      usage:  elmzz=getazz(elem)  

             converts atomic name into atomic number

      example:

            elmzz=getazz("Ga") 

*******************************************************************

   usage:  s = loadstruct(filename)

           reads Wien2k structfile (filename) and returns 
           appropriate structure variable (s)

           filename       is a text variable
           s              is a structure

   Example:

         s = loadstruct("GaN.struct")

******************************************************************* 

   usage:  sr=movealla(s,vec)

        moves all atoms with vector vec

       s       input structure
       vec     vector
       sr      output structure  

  example:

       sr=movealla(s,[0 0 0.5])

*******************************************************************

  usage:  sr=replaceatom(s,ind,atom) 

       replaces an atom ind1 with other atom

       s     input structure
       ind   index of replaced atom
       atom  new atom name 
       sr    output structure 

 example: sr=replaceatom (s,3,"Ge")   

*******************************************************************

  usage:  sr=replaceeqatoms(s,ind,atom,spgn)

          replaces atom ind and all equivalent with other atoms
          
        s        input structure  
        ind      index of removed atoms
        atom     name of new atoms
        spgn     space group name, if not present sym. operations
                 are generated
        sr       output structure

  example: sr=replaceeqatoms (s,3,"Ge","P63mc") 
   
*******************************************************************

   usage:  sr=rmatom(s,ind)

       removes an atom from the structure s and returns all as a structure sr

        sr      output structure
        s       matlab structure used to store input structural information
        ind     a index of a removed atom (position in the list),

  example:
  
         sr=rmatom(s,3)        

*******************************************************************

  usage:  sr=rmeqatoms(s,ind,spgn)

       removes atom ind and all equivalent

       s     input structure 
       ind   index of an atom
       spgn  space group name
       sr    output structure 

  example:  sr=rmeqatoms(s,5,'P63mc')

*******************************************************************

   usage:  savestruct(s,filename,dosym)
          
           saves crystal structure stored in variable s as
           Wien2k structfile.
 
           s              is a structure
           filename       is a text variable
           dosym          1 - with symmetry, 0 - no symmetry

   Example:

         savestruct(s,"GaN.struct",0)

*******************************************************************

   usage:  lea=showequivalent(s,ind,spgn) 

         outputs list of equivalent atoms with atom ind

         s      input structure  
         ind    index of an atom
         spgn   space group name, if not specified sym. operations 
                are degenerated
         lea    output vector containing indexes of eqivalent atoms
     
   example:  el=showequivalent (s2,5,'P63mc')

*******************************************************************

    usage:  showstruct(s)

            displays structure stored in a variable s  
     
  usage:  sr=smultatom(s,ind,spgn)

       creates symmetry equivalent positions for a given atom

       s      input structure 
       ind    index of a multiplied atom 
       spgn   space group symbol, if not specified sym. operations
              are generated 
       sr     output structure

  example:  s2=smultatom (s1,5,"P63mc")

*******************************************************************

   usage:  sr=sshift(s,ind,vec,spgn)

        shifts atom ind with vector vec and all equivalent with
        corresponding equivalents vectors

        s      input structure
        ind    index of an atom  
        vec    displacement vector
        spgn   space group name (optional), if not present
               sym. op. are generated

   example: s2=sshift (s1,5,[ 0.1 0.1 0.1])  

*******************************************************************

   usage:  sr=makeconventional(s)

         converts structure into the conventional form (removes centering)

         sr    output structure
         s     input structure

   example : sr=makeconventional(s)

*******************************************************************

   usage:  sr=makeprimitive(s)

         converts structure to the primitive form (fcc to rhombohedral)

         sr    output structure
         s     input structure

   example : sr=makeprimitive(s)

*******************************************************************

   usage:  sr=makesupercell(s,a)

      creates supercell based on structure s

       s     input structure
       a     new brave lattice in the basis of old lattice vectors
             (row wise), works only with non-centered WIEN structures

   example:  sr=makesupercell(s,[1 0 0; 0 1 0; 0 0 2])

*******************************************************************

   usage:  sr=makesurface(s,n,ind,depth,vac)

      creates surface for a given unitcells

      s       input structure   
      n       normal vector (in lattice coordinates)
      ind     an index of an atom which should be in (0 0 0) 
      depth   thickness of the material  
      vac     thickness of the vacuum layer

   example: sr=makesurface(s,[0 0 1],1,30.0,20.0)

*******************************************************************
