reinitialize
load ./example/cafemol_ligand_hp/cafemol_ligand_hp.movie, object=workobj, format=pdb, multiplex=0

# representation and color
as sphere, (not resn MIY)
alter all, vdw=2.5
rebuild

# do not change atom order in PDB
set retain_order

set ribbon_trace_atoms, 1
set ribbon_radius, 0.7
show_as ribbon, resn MIY 
bg_color white
color blue, chain a
color orange, resn MIY
select exit, (chain a and (i. 20,21,92,95,154,157,158))
color ruby, exit
select entrance, (chain a and (i. 230,231,309-313,366-372))
color red, entrance

# fit frames and smooth the trajectory 
intra_fit (not resn MIY), 1
#       select,  passes, window, first, last, ends
smooth workobj,       3,      5,     1, 2001,   0

# time indicator
set label_position, (0,0,0)
set label_color, black 
set label_size, 36
pseudoatom pa, pos=[12, 45, 0]
hide (pa)

# map states to frames
python
j=1
# total 741 frames
# 101 frames before switching
for i in range(1,1002,10):
  cmd.mset(str(i),str(j))
  state=(i-1)*200
  cmd.mdo(j,"label pa, '"+str(state)+"'")
  j=j+1

# 600 frames after switching and during drug export 
for i in range(1002,1602,1):
  cmd.mset(str(i),str(j))
  state=(i-1)*200
  cmd.mdo(j,"label pa, '"+str(state)+"'")
  j=j+1

# 40 frames after drug export 
for i in range(1611,2002,10):
  cmd.mset(str(i),str(j))
  state=(i-1)*200
  cmd.mdo(j,"label pa, '"+str(state)+"'")
  j=j+1

python end

set_view (\
     0.378969312,   -0.411163270,    0.829050303,\
    -0.784680188,   -0.617686391,    0.052347533,\
     0.490569890,   -0.670377791,   -0.556715667,\
     0.000000000,   -0.000000000, -242.472991943,\
     0.004297256,    0.117418289,    0.011384964,\
   180.192840576,  304.753051758,  -20.000000000 )

frame 1
scene 001, store
mview store, scene=001

frame 701
scene 002, store
mview store, scene=002

frame 702
hide (workobj and resn MIY)
scene 003, store
mview store, scene=003

frame 1

python
for i in range(1,742,1):
  cmd.frame(i)
  cmd.png("t_movgraph/t_"+"%04d"%i, 320, 240, 96, 1)
python end


