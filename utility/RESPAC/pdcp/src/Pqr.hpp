#ifndef RESPAC_PDCP_PQR_H
#define RESPAC_PDCP_PQR_H

#include <string>
#include <vector>

class Pqr
{
    public:

        Pqr(const std::string& fname);
        int size();
        std::string get_atom_name(const int atom_id);
        std::string get_residue_name(const int atom_id);
        int get_residue_id(const int atom_id);
        int get_atom_id(const int atom_id);
        double get_x(const int atom_id);
        double get_y(const int atom_id);
        double get_z(const int atom_id);
        double get_radius(const int atom_id);
        std::vector<double> center_of_mass();
        double get_min_x();
        double get_max_x();
        double get_min_y();
        double get_max_y();
        double get_min_z();
        double get_max_z();
        double get_net_charge();

    private:

        std::vector<double> m_PQR_atom_x;
        std::vector<double> m_PQR_atom_y;
        std::vector<double> m_PQR_atom_z;
        std::vector<double> m_PQR_atom_charge;
        std::vector<double> m_PQR_atom_radius;
        std::vector<std::string> m_PQR_atom_name;
        std::vector<std::string> m_PQR_residue_name;
        std::vector<int> m_PQR_atom_id;
        std::vector<int> m_PQR_residue_id;
        double m_net_charge;
        int m_n_atom;
};

#endif //RESPAC_PDCP_PQR_H
