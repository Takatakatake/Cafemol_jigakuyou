#!/bin/bash

INPUTFILE=./example/cafemol_ligand_hp/cafemol_ligand_hp.movie
OUTPUTFILE=./cafe_movie.avi
SCRIPT_PATH=./utility
TEMP_PATH=./tmp
DEL_TEMP_FILES=0

mk_temp_path=0
if [ ! -r $TEMP_PATH ]; then 
   mkdir $TEMP_PATH
   mk_temp_path=1
fi
rm -rf $TEMP_PATH/t_movgraph
mkdir $TEMP_PATH/t_movgraph


# Needs PyMol 1.2 (or above) installed
# See http://www.pymol.org
sed -e "s#^load .*, object#load $INPUTFILE, object#" -e "s#t_movgraph#$TEMP_PATH/t_movgraph#" \
  $SCRIPT_PATH/cafe_movie_making_gengraph.pml \
  > $TEMP_PATH/cafe_movie_making_gengraph_myjob.pml 
pymol -c -u $TEMP_PATH/cafe_movie_making_gengraph_myjob.pml 



# Needs MPlayer package installed
# See http://www.mplayerhq.hu
mencoder -mc 0 -noskip -skiplimit 0 -ovc lavc -lavcopts \
vcodec=mpeg4:vhq:trell:mbd=2:vmax_b_frames=1:v4mv:vb_strategy=0:o=luma_elim_threshold=0:\
o=chroma_elim_threshold=0:cmp=6:subcmp=6:precmp=6:predia=3:dia=3:vme=4:vqscale=1 \
mf://$TEMP_PATH/t_movgraph/*.png -mf type=png:fps=24 -o $OUTPUTFILE

if [ $DEL_TEMP_FILES -eq 1 ]; then 
   rm -rf $TEMP_PATH/t_movgraph $TEMP_PATH/cafe_movie_making_gengraph_myjob.pml
   if [ $mk_temp_path -eq 1 ]; then rm -rf $TEMP_PATH; fi
fi
