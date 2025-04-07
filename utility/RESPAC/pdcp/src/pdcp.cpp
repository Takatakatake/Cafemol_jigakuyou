#include <iostream>
#include "Input.hpp"
#include "Parameter.hpp"
#include "Pqr.hpp"
#include "Pdc.hpp"
#include "Dx.hpp"
#include "util.hpp"
#include "Charge.hpp"

int main(int argc, char** argv)
{
    std::cerr << std::endl;
    std::cerr << "\t*******************************************" << std::endl;
    std::cerr << "\t*                 ePdcp 0.1               *" << std::endl;
    std::cerr << "\t* originally written by Tsuyoshi Terakawa *" << std::endl;
    std::cerr << "\t*           modified by Toru Niina        *" << std::endl;
    std::cerr << "\t*******************************************" << std::endl;
    std::cerr << std::endl;
 
    Input input(argc, argv);
    Parameter parameter(input.get_ifname());
    Pqr pqr(input.get_pqr_fname());
    Pdc pdc(pqr, input);
    Dx dx_pot(input.get_pot_fname());
    Dx dx_vol(input.get_vol_fname());
 
    dx_pot = dx_pot * dx_vol;
    double net_charge = pqr.get_net_charge();

    std::cerr << "Net Charge: " << net_charge << std::endl;
 
    if(input.get_kai_fname() == "")
    {
        calculate_pdc(dx_pot, pdc, parameter, input.get_ofname(), net_charge);
    }
    else
    {
        //cout << input.get_kai_fname() << endl;
        Charge charge(input.get_kai_fname());
        calculate_chi(dx_pot, pdc, charge, parameter);
    }
 
    return 0;
}
