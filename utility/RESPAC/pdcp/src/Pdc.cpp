#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sys/stat.h>
#include <boost/lexical_cast.hpp>
#include "Pdc.hpp"
using namespace boost;

Pdc::Pdc(Pqr& pqr, Input& input)
{
    // Make map for judge
    init_map();
 
    if(input.get_res_fname() != "")
    {
        init_specified_residue(input.get_res_fname());
    }
 
    m_size = 0;
    for(int iatom = 0; iatom < pqr.size(); iatom++)
    {
        std::string atom_name    = pqr.get_atom_name(iatom);
        std::string residue_name = pqr.get_residue_name(iatom);
        int residue_id = pqr.get_residue_id(iatom);
        
        if((input.get_site() == "All" && judge_charge_particle(atom_name)) || 
           (input.get_site() == "Charge" && judge_charge_particle(atom_name, residue_name)))
        {
            if(input.get_res_fname() == "" || judge_charge_particle(residue_id))
            {
                m_size++;
                m_residue_id.push_back(residue_id);
                m_pos_x.push_back(pqr.get_x(iatom));
                m_pos_y.push_back(pqr.get_y(iatom));
                m_pos_z.push_back(pqr.get_z(iatom));
            }
        } 
    }
    std::cerr << "The # of charge sites: " << m_size << std::endl;
    std::cerr << std::endl;
}

void Pdc::init_map()
{
    m_charged_residue["ARG"] = true;
    m_charged_residue["LYS"] = true;
    m_charged_residue["HIS"] = true;
    m_charged_residue["ASP"] = true;
    m_charged_residue["GLU"] = true;
    return;
}

void Pdc::init_specified_residue(const std::string& ifname)
{
    struct stat st;
    int ret = stat(ifname.c_str(), &st);
    if(ret != 0)
    {
        std::cerr << "The file \"" << ifname << "\" does not exist: File_name = "
                  << __FILE__ << " Line = " << __LINE__  << std::endl;
        exit(1);
    }
    
    std::ifstream ifs(ifname.c_str());
    std::string line;
    while(std::getline(ifs, line))
    {
        m_specified_residue[lexical_cast<int>(line)] = true;
    }
    ifs.close();
    return;
}

bool Pdc::judge_charge_particle(const std::string& atom_name)
{
    return (atom_name == "CA");
}

bool Pdc::judge_charge_particle(const std::string& atom_name,
                                const std::string& residue_name)
{
    return (m_charged_residue.find(residue_name) != m_charged_residue.end() && 
            atom_name == "CA");
}

bool Pdc::judge_charge_particle(const int residue_id)
{
    return (m_specified_residue.find(residue_id) != m_specified_residue.end());
}
