default:

install:
	cd downloads; tar xvf x3dna*.tar.gz
	cd downloads; unzip USER-3SPN2.zip
	cd downloads/USER-3SPN2/DSIM_ICNF/; make all
	cd downloads/USER-3SPN2/utils; chmod +x *
run:
	cd utils; make run
clean:
	rm -f *.pdb strand*.ninfo
	cd utils; make clean
purge:  clean
	rm -rf downloads/USER-3SPN2/ downloads/x3dna*/
