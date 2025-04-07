#include <iostream>
#include <string>
#include <cstdlib>
#include <boost/program_options.hpp>
#include "Input.hpp"
using namespace std;

namespace po = boost::program_options;

Input::Input(int ac, char** av) {
  po::options_description desc("Allowed options");

  desc.add_options()
    ("help", "produce help message")
    ("pqr", po::value<string>(&m_pqr_fname)->default_value("protein.pqr"), "PQR file name")
    ("ofname", po::value<string>(&m_ofname)->default_value("output"), "Output file name")
    ("dbox", po::value<double>(&m_dbox)->default_value(2.0), "Minimum distance between protein and box")
    ("dx", po::value<double>(&m_dx)->default_value(1.0), "Grid resolution")
    ("r_probe", po::value<double>(&m_r_probe)->default_value(1.4), "Radius of probe")
    ;

  po::variables_map vm;
  po::store(po::parse_command_line(ac, av, desc), vm);
  po::notify(vm);

  if (vm.count("help")) {
    cerr << desc << "\n";
    exit(1);
  }

  cerr << "PQR file name: " << m_pqr_fname << endl;
  cerr << "Output file name: " << m_ofname << endl;
  cerr << "Minimum distance between protein and box: " << m_dbox << endl;
  cerr << "Grid resolution: " << m_dx << endl;
  cerr << "Radius of probe: " << m_r_probe << endl << endl;
  
}

string Input::get_ofname() {
  return m_ofname;
}

string Input::get_pqr_fname() {
  return m_pqr_fname;
}

double Input::get_dbox() {
  return m_dbox;
}

double Input::get_dx() {
  return m_dx;
}

double Input::get_r_probe() {
  return m_r_probe;
}
