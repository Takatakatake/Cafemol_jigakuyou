.SUFFIXES : .o .F90
TARGET = DNAbackmap

OBJ   = const_para.o ifmod.o util_bestfit.o output_pdb.o in_out_file.o\
	end_modeling.o main_modeling.o check_bp.o read_pdb.o \
	read_references.o DNAbackmap.o
MOD   = const_para.mod interface_modeling.mod

FC = gfortran
FC_UTIL = gfortran
#OPT = -xHOST -O3 -no-prec-div -static
OPT = 
#LIB = 


.F90.o:
	${FC} -c $<

${TARGET} : ${OBJ}
	${FC} -o $@ ${OBJ} ${OPT} 

.PHONY: clean
clean:
	rm -f ${TARGET} ${OBJ} ${MOD}
