#include <iostream>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <sys/stat.h>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/lexical_cast.hpp>
#include "Pqr.hpp"
using namespace boost;
using namespace boost::algorithm;

Pqr::Pqr(const std::string& fname)
{
    struct stat st;
    int ret = stat(fname.c_str(), &st);
    if(ret != 0)
    {
        std::cerr << "The file \"" << fname << "\" does not exist: File_name = "
                  << __FILE__ << " Line = " << __LINE__  << std::endl;
        exit(1);
    }
 
    std::ifstream ifs(fname.c_str());

    std::string line;
    std::vector<std::string> columns;
    while(std::getline(ifs, line))
    {
        //split(columns, line, is_any_of(" "), token_compress_on);
        if(line.size() < 4) continue;
        if(line.substr(0, 4) == "ATOM")
        {
            m_PQR_atom_x.push_back(lexical_cast<double>(trim_copy(line.substr(30, 8))));
            m_PQR_atom_y.push_back(lexical_cast<double>(trim_copy(line.substr(38, 8))));
            m_PQR_atom_z.push_back(lexical_cast<double>(trim_copy(line.substr(46, 8))));
            m_PQR_atom_charge.push_back(lexical_cast<double>(trim_copy(line.substr(54, 8))));
            m_PQR_atom_radius.push_back(lexical_cast<double>(trim_copy(line.substr(62, 8))));
            m_PQR_atom_name.push_back(trim_copy(line.substr(12, 4)));
            m_PQR_residue_name.push_back(trim_copy(line.substr(17, 4)));
            m_PQR_atom_id.push_back(lexical_cast<int>(trim_copy(line.substr(5, 6))));
            m_PQR_residue_id.push_back(lexical_cast<int>(trim_copy(line.substr(20, 6))));
        }
    }
 
    m_n_atom = m_PQR_atom_x.size();
 
    m_net_charge = 0.0;
    for(std::size_t i = 0; i < m_PQR_atom_charge.size(); i++)
    {
        m_net_charge += m_PQR_atom_charge[i];
    }
 
    //for (vector<int>::iterator it = m_PQR_residue_id.begin(); it != m_PQR_residue_id.end(); ++it) {
    //  cout << *it << endl;
    //}
}

int Pqr::size()
{
    return m_n_atom;
}

std::string Pqr::get_atom_name(const int atom_id)
{
    return m_PQR_atom_name[atom_id];
}

std::string Pqr::get_residue_name(const int atom_id)
{
    return m_PQR_residue_name[atom_id];
}

int Pqr::get_residue_id(const int atom_id)
{
    return m_PQR_residue_id[atom_id];
}

int Pqr::get_atom_id(const int atom_id)
{
    return m_PQR_atom_id[atom_id];
}

double Pqr::get_x(const int atom_id)
{
    return m_PQR_atom_x[atom_id];
}

double Pqr::get_y(const int atom_id)
{
    return m_PQR_atom_y[atom_id];
}

double Pqr::get_z(const int atom_id)
{
    return m_PQR_atom_z[atom_id];
}

double Pqr::get_radius(const int atom_id)
{
    return m_PQR_atom_radius[atom_id];
}

std::vector<double> Pqr::center_of_mass()
{
    std::vector<double> center(3, 0.0);
    int n = 0;
 
    for(std::size_t i = 0; i < m_PQR_atom_x.size(); ++i)
    {
        center[0] += m_PQR_atom_x[i];
        center[1] += m_PQR_atom_y[i];
        center[2] += m_PQR_atom_z[i];
        n++;
    }
 
    for(int i = 0; i < 3; ++i)
    {
        center[i] /= static_cast<double>(n);
    }

    return center;
}

double Pqr::get_min_x()
{
    double coord = m_PQR_atom_x[0];
 
    for(std::size_t i = 1; i < m_PQR_atom_x.size(); ++i)
    {
        if(coord > m_PQR_atom_x[i])
        {
            coord = m_PQR_atom_x[i];
        }
    }
 
    return coord;
}

double Pqr::get_max_x()
{
    double coord = m_PQR_atom_x[0];
 
    for(std::size_t i = 1; i < m_PQR_atom_x.size(); ++i)
    {
        if(coord < m_PQR_atom_x[i])
        {
            coord = m_PQR_atom_x[i];
        }
    }
 
    return coord;
}

double Pqr::get_min_y()
{
    double coord = m_PQR_atom_y[0];
 
    for(std::size_t i = 1; i < m_PQR_atom_y.size(); ++i)
    {
        if(coord > m_PQR_atom_y[i])
        {
            coord = m_PQR_atom_y[i];
        }
    }
 
    return coord;
}

double Pqr::get_max_y()
{
    double coord = m_PQR_atom_y[0];
 
    for(std::size_t i = 1; i < m_PQR_atom_y.size(); ++i)
    {
        if(coord < m_PQR_atom_y[i])
        {
            coord = m_PQR_atom_y[i];
        }
    }
 
    return coord;
}

double Pqr::get_min_z()
{
    double coord = m_PQR_atom_z[0];
 
    for(std::size_t i = 1; i < m_PQR_atom_z.size(); ++i)
    {
        if(coord > m_PQR_atom_z[i])
        {
            coord = m_PQR_atom_z[i];
        }
    }
 
    return coord;
}

double Pqr::get_max_z()
{
    double coord = m_PQR_atom_z[0];
 
    for(std::size_t i = 1; i < m_PQR_atom_z.size(); ++i)
    {
        if(coord < m_PQR_atom_z[i])
        {
            coord = m_PQR_atom_z[i];
        }
    }
 
    return coord;
}

double Pqr::get_net_charge()
{
    return m_net_charge;
}
