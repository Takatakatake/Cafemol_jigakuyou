#ifndef _GRID_H
#define _GRID_H

#include "Pqr.hpp"
#include "Input.hpp"

class Grid {
 private:
  vector<double> m_coord_x;
  vector<double> m_coord_y;
  vector<double> m_coord_z;
  int m_nx;
  int m_ny;
  int m_nz;
 public:
  Grid(Pqr& pqr, Input& input);
  int get_nx();
  int get_ny();
  int get_nz();
  int size();
  double get_x(const int& id);
  double get_y(const int& id);
  double get_z(const int& id);
};

#endif // #ifndef _GRID_H
