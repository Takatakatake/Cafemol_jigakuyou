#ifndef RESPAC_PDCP_PDC_H
#define RESPAC_PDCP_PDC_H

#include <vector>
#include <map>
#include <string>
#include "Pqr.hpp"
#include "Input.hpp"

class Pdc
{
    public:

        Pdc(Pqr& pqr, Input& input);

        int size();
        double get_x(const int id);
        double get_y(const int id);
        double get_z(const int id);
        int get_residue_id(const int id);


    private:

        std::vector<double> m_pos_x;
        std::vector<double> m_pos_y;
        std::vector<double> m_pos_z;
        std::vector<int> m_residue_id;
        std::map<std::string, bool> m_charged_residue;
        std::map<int, bool> m_specified_residue;
        int m_size;

        void init_map();
        void init_specified_residue(const std::string& ifname);
        bool judge_charge_particle(const std::string& atom_name);
        bool judge_charge_particle(const std::string& atom_name,
                                   const std::string& residue_name);
        bool judge_charge_particle(const int residue_id);
};

inline int Pdc::size()
{
    return m_size;
}

inline double Pdc::get_x(const int id)
{
    return m_pos_x[id];
}

inline double Pdc::get_y(const int id)
{
    return m_pos_y[id];
}

inline double Pdc::get_z(const int id)
{
    return m_pos_z[id];
}

inline int Pdc::get_residue_id(const int id)
{
    return m_residue_id[id];
}

#endif //RESPAC_PDCP_PDC_H
