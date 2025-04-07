#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include "Charge.hpp"
using namespace boost;
using namespace boost::algorithm;

Charge::Charge(const std::string& fname)
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

    std::vector<std::string> fnames;
    std::string line;
    while(std::getline(ifs, line))
    {
        fnames.push_back(line);
    }

    ifs.close();

    for(std::vector<std::string>::iterator it = fnames.begin();
        it != fnames.end(); ++it)
    {
        ret = stat((*it).c_str(), &st);
        if(ret != 0)
        {
            std::cerr << "The file \"" << *it << "\" does not exist: File_name = "
                      << __FILE__ << " Line = " << __LINE__  << std::endl;
            exit(1);
        }

        std::ifstream ifs_local((*it).c_str());
        std::vector<double> tmp_v;
        std::vector<std::string> columns;
        while(std::getline(ifs_local, line))
        {
            split(columns, line, is_any_of(" "), token_compress_on);
            tmp_v.push_back(lexical_cast<double>(columns[1]));
        }

        m_charges.push_back(tmp_v);

        ifs_local.close();
    }

    m_n_charge    = m_charges[0].size();
    m_n_structure = m_charges.size();
}

double Charge::operator()(const int structure_id, const int charge_id)
{
    return m_charges[structure_id][charge_id];
}

int Charge::get_n_charge()
{
    return m_n_charge;
}

int Charge::get_n_structure()
{
    return m_n_structure;
}
