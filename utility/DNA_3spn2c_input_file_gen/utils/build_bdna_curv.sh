#!/bin/sh

# This tool is originally built by de Pablo's group and modified to generate
# input files for CafeMol.

if [ $# -ne 2 ]; then
    echo "Usage: $0 <sequence file> <strand info file>"
    exit 1
fi

ICNF=../downloads/USER-3SPN2/DSIM_ICNF/icnf.exe
UTILS=../downloads/USER-3SPN2/utils

echo "Making parameter file"
python2 $UTILS/make_bp_params.py $1

echo "Running X3DNA"
# export PATH=$PATH:../downloads/x3dna-v2.2/bin
export PATH=$PATH:../downloads/x3dna-v2.3/bin
# export X3DNA=../downloads/x3dna-v2.2
export X3DNA=../downloads/x3dna-v2.3
x3dna_utils cp_std BDNA
rebuild -atomic bp_step.par atomistic.pdb
mkdir basis
mv Atomic* basis/

echo "Mapping to CG coordinates"
python2 $UTILS/pdb2cg_dna.py atomistic.pdb
mv in00_conf.xyz bdna_curv.xyz

# echo "Building LAMMPS input file"
${ICNF} $1 1 1 . 0 1> /dev/null

echo "Replacing atoms in configuration file"
$UTILS/replace_atoms.sh conf_lammps.in dna_conf.in bdna_curv_conf.in

echo "Making list files"
python2 $UTILS/make_list_files.py bdna_curv_conf.in

echo "Making cafemol files"
./gen_cafemol_files.py $2
