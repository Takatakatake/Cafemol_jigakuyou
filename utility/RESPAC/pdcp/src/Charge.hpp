#ifndef RESPAC_PDCP_CHARGE_H
#define RESPAC_PDCP_CHARGE_H

#include <string>
#include <vector>

class Charge
{
    public:

        Charge(const std::string& fname);

        double operator()(const int structure_id, const int charge_id);
        int get_n_charge();
        int get_n_structure();

    private:

        std::vector<std::vector<double> > m_charges;
        int m_n_charge;
        int m_n_structure;
};

#endif //RESPAC_PDCP_CHARGE_H
