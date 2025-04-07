#include <iostream>
#include "Grid.hpp"
using namespace std;

Grid::Grid(Pqr& pqr, Input& input) {
  
  double max_x = pqr.get_max_x() + input.get_dbox();
  double min_x = pqr.get_min_x() - input.get_dbox();
  double max_y = pqr.get_max_y() + input.get_dbox();
  double min_y = pqr.get_min_y() - input.get_dbox();
  double max_z = pqr.get_max_z() + input.get_dbox();
  double min_z = pqr.get_min_z() - input.get_dbox();
  double dx = input.get_dx();

  double x = min_x;
  m_nx = 0;
  while (1) {
    m_coord_x.push_back(x);
    m_nx++;
    if (x >= max_x) {
      break;
    }
    x += dx;
  }
  
  double y = min_y;
  m_ny = 0;
  while (1) {
    m_coord_y.push_back(y);
    m_ny++;
    if (y >= max_y) {
      break;
    }
    y += dx;
  }

  double z = min_z;
  m_nz = 0;
  while (1) {
    m_coord_z.push_back(z);
    m_nz++;
    if (z >= max_z) {
      break;
    }
    z += dx;
  }
  
}

int Grid::get_nx () {
  return m_nx;
}

int Grid::get_ny () {
  return m_ny;
}

int Grid::get_nz () {
  return m_nz;
}

int Grid::size () {
  return m_nx * m_ny * m_nz;
}

double Grid::get_x (const int& id) {
  return m_coord_x[id];
}

double Grid::get_y (const int& id) {
  return m_coord_y[id];
}

double Grid::get_z (const int& id) {
  return m_coord_z[id];
}
