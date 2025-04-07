#ifndef _INPUT_H
#define _INPUT_H

#include <string>
using namespace std;

class Input {
 private:
  string m_pqr_fname;
  string m_ofname;
  double m_dbox;
  double m_dx;
  double m_r_probe;
 public:
  Input(int ac, char** av);
  string get_pqr_fname();
  string get_ofname();
  double get_dbox();
  double get_dx();
  double get_r_probe();
};

#endif // #ifndef _INPUT_H
