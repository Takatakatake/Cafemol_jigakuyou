#include <iostream>
#include <string>
#include <cstdlib>
#include <boost/program_options.hpp>
#include "Input.hpp"

namespace po = boost::program_options;

Input::Input(int ac, char** av)
{
    po::options_description desc("Allowed options");
 
    desc.add_options()
        ("help", "produce help message")
        ("ifname",
         po::value<std::string>(&m_ifname)->default_value("input"), 
         "Input file name")
        ("ofname",
         po::value<std::string>(&m_ofname)->default_value("output"),
         "Output file name")
        ("pqr",
         po::value<std::string>(&m_pqr_fname)->default_value("protein.pqr"),
         "PQR file name")
        ("pot",
         po::value<std::string>(&m_pot_fname)->default_value("potential.dx"),
         "Potential file name")
        ("vol",
         po::value<std::string>(&m_vol_fname)->default_value("volume.dx"),
         "Volume file name")
        ("kai",
         po::value<std::string>(&m_kai_fname)->default_value(""),
         "File name for kai calculation")
        ("site",
         po::value<std::string>(&m_site)->default_value("All"),
         "The site where charges are put [All or Charge]")
        ("residue",
         po::value<std::string>(&m_res_fname)->default_value(""),
         "Residue file name");
 
    po::variables_map vm;
    po::store(po::parse_command_line(ac, av, desc), vm);
    po::notify(vm);
 
    if(vm.count("help"))
    {
        std::cerr << desc << std::endl;
        std::exit(1);
    }
 
    if(m_site != "All" && m_site != "Charge")
    {
        std::cerr << "\"site\" should be \"All\" or \"Charge\"" << std::endl;
        throw std::invalid_argument("invalid site value");
    }
 
    std::cerr << "Info >> Input" << std::endl;
    std::cerr << "Input file name       : " << m_ifname    << std::endl;
    std::cerr << "Output file name      : " << m_ofname    << std::endl;
    std::cerr << "PQR file name         : " << m_pqr_fname << std::endl;
    std::cerr << "Potential file name   : " << m_pot_fname << std::endl;
    std::cerr << "Volume file name      : " << m_vol_fname << std::endl;

    if(!m_res_fname.empty())
    {
        std::cerr << "Residue file name     : " << m_res_fname << std::endl;
    }

    std::cerr << "The site where charges are put : " << m_site << std::endl;
    std::cerr << std::endl;
}
