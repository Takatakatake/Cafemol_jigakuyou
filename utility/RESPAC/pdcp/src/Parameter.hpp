#ifndef RESPAC_PDCP_PARAMETER_H
#define RESPAC_PDCP_PARAMETER_H

#include <string>

class Parameter
{
    public:

        Parameter(const std::string& fname);

        double get_temp()    const {return m_temp;}
        double get_debye()   const {return m_debye;}
        double get_diel()    const {return m_diel;}
        double get_penalty() const {return m_penalty;}
        double get_net_charge_constraint() const
                                   {return m_net_charge_constraint;}

    private:

        int m_cpu;
        int m_points;
        double m_diel;
        double m_temp;
        double m_debye;
        double m_penalty;
        double m_net_charge_constraint;
};

#endif //RESPAC_PDCP_PARAMETER_H
