#!/usr/bin/env python

import sys, glob, shutil
from os import path

if len(sys.argv) != 2:
   print ('Usage: SCRIPT [dir name]')
   sys.exit(2)

dirname = sys.argv[1]
if not dirname[:1] in ('.','/'):
   dirname = './' + dirname

if (not path.exists(dirname) or not path.isdir(dirname)):
   print ('directory %s does not exist' % dirname)
   sys.exit(2)

exts = ('ts','pdb','dcd','psf','crd','velo','movie',
        'data','rst','vdcd','ninfo','opt','rep')

for ex in exts:
   pattern = '*.%s' % ex
   for file in glob.glob(path.join(dirname,pattern)):
      newfile = file + '_example'
      shutil.move(file, newfile)
   
