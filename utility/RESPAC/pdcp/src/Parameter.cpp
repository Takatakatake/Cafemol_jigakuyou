#include <iostream>
#include <fstream>
#include <vector>
#include <stdexcept>
#include <cstdlib>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include "Parameter.hpp"

Parameter::Parameter(const std::string& fname)
{
    std::ifstream ifs(fname);
    if(!ifs.good())
    {
        std::cerr << "The file \"" << fname << "\" cannot open"<< std::endl;
        throw std::invalid_argument("file open error");
    }
 
    std::string line;
    std::vector<std::string> columns;
    while(std::getline(ifs, line))
    {
        boost::algorithm::split(
                    columns,
                    line,
                    boost::is_any_of(" "),
                    boost::algorithm::token_compress_on
                );

        if(columns.size() != 2)
        {
            std::cerr << "Column size must be 2. " << std::endl;
            std::cerr << "line : " << line << std::endl;
            throw std::invalid_argument("invalid file format");
        }

        if(columns[0] == "diel")
        {
            m_diel = std::stod(columns[1]);
        }
        else if(columns[0] == "temp")
        {
            m_temp = std::stod(columns[1]);
        }
        else if(columns[0] == "debye")
        {
            m_debye = std::stod(columns[1]);
        }
        else if(columns[0] == "penalty")
        {
            m_penalty = std::stod(columns[1]);
        }
        else if(columns[0] == "points")
        {
            m_points = std::stoi(columns[1]);
        }
        else if(columns[0] == "CPU")
        {
            m_cpu = std::stoi(columns[1]);
        }
        else if(columns[0] == "NetCharge")
        {
            m_net_charge_constraint = std::stod(columns[1]);
        }
        else
        {
            std::cerr << "Unknown line : " << line << std::endl;
        }
    }
    ifs.close();
 
    std::cerr << "Info >> Parameters" << std::endl;
    std::cerr << "Dielectric constant   : " << m_diel
              << std::endl;
    std::cerr << "Temperature           : " << m_temp
              << std::endl;
    std::cerr << "Debye length          : " << m_debye
              << std::endl;
    std::cerr << "Penalty               : " << m_penalty
              << std::endl;
    std::cerr << "Net charge constraint : " << m_net_charge_constraint
              << std::endl;
    std::cerr << "Points                : " << m_points
              << std::endl;
    std::cerr << "CPU                   : " << m_cpu
              << std::endl;
    std::cerr << std::endl;

}
