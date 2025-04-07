#ifndef RESPAC_PDCP_UTIL_H
#define RESPAC_PDCP_UTIL_H

#include <string>
#include "Dx.hpp"
#include "Pdc.hpp"
#include "Parameter.hpp"
#include "Charge.hpp"

void calculate_pdc(Dx& dx_pot,
                   Pdc& pdc,
                   Parameter& param,
                   const std::string& ofname,
                   double net_charge);

void calculate_chi(Dx& dx_pot,
                   Pdc& pdc,
                   Charge& charge,
                   Parameter& param);

#endif //RESPAC_PDCP_UTIL_H
