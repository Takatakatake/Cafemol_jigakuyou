#include <iostream>
#include <fstream>
#include <map>
#include <cmath>
#include <limits.h>
#include "util.hpp"
using namespace std;

void decide_surface(Pqr& pqr, Grid& grid, Input& input) {
  vector<int> surface_flag( grid.size(), 1 );

  int ipoint = 0;
  for (int ix = 0; ix < grid.get_nx(); ++ix) {
    for (int iy = 0; iy < grid.get_ny(); ++iy) {
      for (int iz = 0; iz < grid.get_nz(); ++iz) {
	double grid_x = grid.get_x(ix);
	double grid_y = grid.get_y(iy);
	double grid_z = grid.get_z(iz);

	for (int iatom = 0; iatom < pqr.size(); ++iatom) {
	  double pqr_x = pqr.get_x(iatom);
	  double pqr_y = pqr.get_y(iatom);
	  double pqr_z = pqr.get_z(iatom);
	  
	  double dist = sqrt( (grid_x - pqr_x) * (grid_x - pqr_x) + 
			      (grid_y - pqr_y) * (grid_y - pqr_y) + 
			      (grid_z - pqr_z) * (grid_z - pqr_z) );
	  
	  if (pqr.get_radius(iatom) + input.get_r_probe() >= dist) { 
	    surface_flag[ipoint] = 0;
	  }
	  
	}
	
	ipoint++;
	
      }
    }
  }

  ipoint = 0;
  map<int, int> residues;
  for (int ix = 0; ix < grid.get_nx(); ++ix) {
    for (int iy = 0; iy < grid.get_ny(); ++iy) {
      for (int iz = 0; iz < grid.get_nz(); ++iz) {

	if (surface_flag[ipoint]) {
	  double min_atom = -1;
	  double min_dist = INT_MAX;
	
	  double grid_x = grid.get_x(ix);
	  double grid_y = grid.get_y(iy);
	  double grid_z = grid.get_z(iz);
	  
	  for (int iatom = 0; iatom < pqr.size(); ++iatom) {
	    double pqr_x = pqr.get_x(iatom);
	    double pqr_y = pqr.get_y(iatom);
	    double pqr_z = pqr.get_z(iatom);
	    
	    double dist = sqrt( (grid_x - pqr_x) * (grid_x - pqr_x) + 
				(grid_y - pqr_y) * (grid_y - pqr_y) + 
				(grid_z - pqr_z) * (grid_z - pqr_z) );
	    
	    if (min_dist > dist) {
	      min_dist = dist;
	      min_atom = iatom;
	    }
	    
	  }

	  residues[pqr.get_residue_id(min_atom)] = 1;
	  //cout << pqr.get_residue_id(min_atom) << " " << min_dist << endl;
	  
	} // if
	
	ipoint++;
      } // iz
    } // iy
  } // ix

  ofstream ofs(input.get_ofname().c_str());
  
  map<int, int>::iterator it = residues.begin();
  while( it != residues.end() ) {
    ofs  << (*it).first << endl;
    cerr << (*it).first << ",";
    ++it;
  }

  cerr << endl;
  
}
