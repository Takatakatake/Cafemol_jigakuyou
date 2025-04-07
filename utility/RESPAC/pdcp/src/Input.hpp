#ifndef RESPAC_PDCP_INPUT_H
#define RESPAC_PDCP_INPUT_H

#include <string>

class Input
{
    public:

        Input(int ac, char** av);
        const std::string get_ifname()    const {return m_ifname;}
        const std::string get_ofname()    const {return m_ofname;}
        const std::string get_pqr_fname() const {return m_pqr_fname;}
        const std::string get_site()      const {return m_site;}
        const std::string get_pot_fname() const {return m_pot_fname;}
        const std::string get_vol_fname() const {return m_vol_fname;}
        const std::string get_kai_fname() const {return m_kai_fname;}
        const std::string get_res_fname() const {return m_res_fname;}

    private:

        std::string m_ifname;
        std::string m_ofname;
        std::string m_pqr_fname;
        std::string m_pot_fname;
        std::string m_vol_fname;
        std::string m_kai_fname;
        std::string m_res_fname;
        std::string m_site;
};

#endif //RESPAC_PDCP_INPUT_H
