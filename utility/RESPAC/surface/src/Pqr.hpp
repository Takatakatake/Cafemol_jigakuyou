#ifndef _PQR_H
#define _PQR_H

#include <string>
#include <vector>
using namespace std;

class Pqr {
 private:
  vector<double> m_PQR_atom_x;
  vector<double> m_PQR_atom_y;
  vector<double> m_PQR_atom_z;
  vector<double> m_PQR_atom_charge;
  vector<double> m_PQR_atom_radius;
  vector<string> m_PQR_atom_name;
  vector<string> m_PQR_residue_name;
  vector<int> m_PQR_atom_id;
  vector<int> m_PQR_residue_id;
  int m_n_atom;

 public:
  Pqr(const string& fname);
  int size();
  string get_atom_name(const int& atom_id);
  string get_residue_name(const int& atom_id);
  int get_residue_id(const int& atom_id);
  int get_atom_id(const int& atom_id);
  double get_x(const int& atom_id);
  double get_y(const int& atom_id);
  double get_z(const int& atom_id);
  double get_radius(const int& atom_id);
  vector<double> center_of_mass();
  double get_min_x();
  double get_max_x();
  double get_min_y();
  double get_max_y();
  double get_min_z();
  double get_max_z();
};

#endif // #ifndef _PQR_H
