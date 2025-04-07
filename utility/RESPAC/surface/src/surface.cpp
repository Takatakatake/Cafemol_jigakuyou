#include <iostream>
#include "Input.hpp"
#include "Pqr.hpp"
#include "Grid.hpp"
#include "util.hpp"
using namespace std;

int main (int argc, char** argv) {

  cerr << endl;
  cerr << "\t******************************" << endl;
  cerr << "\t*        surface 0.0         *" << endl;
  cerr << "\t******************************" << endl << endl;

  Input input(argc, argv);
  Pqr pqr( input.get_pqr_fname() );
  Grid grid(pqr, input);

  decide_surface(pqr, grid, input);
  
  return 0;
}
