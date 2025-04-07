#include <iostream>
#include <fstream>
#include <cmath>
/////#include <gmd.h>
/////#include <lavd.h>
/////#include <laslv.h>

#ifdef PDCP_USE_FMATH
#include "fmath/fmath.hpp"
#endif

#include "Eigen/Dense"
#include "util.hpp"

void calculate_pdc(Dx& dx_pot,
                   Pdc& pdc,
                   Parameter& param,
                   const std::string& ofname,
                   double net_charge)
{
    // Make table
    int i_phi  = 0;
    int i_table(0);
    std::vector<double> phi_table(dx_pot.get_n_data(), 0e0);
    std::vector<int>    x_table(dx_pot.get_n_data(), 0);
    std::vector<int>    y_table(dx_pot.get_n_data(), 0);
    std::vector<int>    z_table(dx_pot.get_n_data(), 0);
    for(int ix = 0; ix < dx_pot.get_grid_x(); ++ix)
    {
        for(int iy = 0; iy < dx_pot.get_grid_y(); ++iy)
        {
            for(int iz = 0; iz < dx_pot.get_grid_z(); ++iz)
            {
                double data = dx_pot.get_data(i_phi);
                if(std::abs(data) > 10e-16)
                {
                    phi_table[i_table] = data;
                    x_table[i_table] = ix;
                    y_table[i_table] = iy;
                    z_table[i_table] = iz;
                    ++i_table;
                }
                i_phi++;
            }
        }
    }

    int n_eff = pdc.size();
    int n_grid_point = i_table;

    std::cerr << "The # of grit point used to fit: "
              << n_grid_point << std::endl;
    std::cerr << std::endl;

    // Initialize matrix (fill zero)
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n_eff, n_eff);
    Eigen::VectorXd b = Eigen::VectorXd::Zero(n_eff);

    // Filling matrix A
    double tmp_val;
    // double factor  = 1.87944 * param.get_temp();
    double factor  = 559.219;
    double invd    = 1.0 / param.get_debye();
    double diel    = param.get_diel();
    double penalty = param.get_penalty();
    double net_charge_constraint = param.get_net_charge_constraint();

    std::cout << "factor: " << factor << std::endl;
    std::cout << "invd: " << invd << std::endl;
    std::cout << "diel: " << diel << std::endl;
    std::cout << "penalty: " << penalty << std::endl;
    std::cout << "net_charge_constraint: " << net_charge_constraint << std::endl;

    double factor2(factor * factor);
    double diel2(diel * diel);
    double temp_coef(factor2 / diel2);

    for(int i = 0; i < n_eff; i++)
    {
        double x_chr_i = pdc.get_x(i); 
        double y_chr_i = pdc.get_y(i); 
        double z_chr_i = pdc.get_z(i); 

        // i == j case
        for(int k(0); k < n_grid_point; ++k)
        {
            double x_pot = dx_pot.get_x(x_table[k]);
            double y_pot = dx_pot.get_y(y_table[k]);
            double z_pot = dx_pot.get_z(z_table[k]);
            double disi2 = (x_pot - x_chr_i)*(x_pot - x_chr_i) +
                           (y_pot - y_chr_i)*(y_pot - y_chr_i) +
                           (z_pot - z_chr_i)*(z_pot - z_chr_i);
#ifdef PDCP_USE_FMATH
            A(i, i) += fmath::expd(-2e0 * invd * std::sqrt(disi2)) / disi2;
#else
            A(i, i) += std::exp(-2e0 * invd * std::sqrt(disi2)) / disi2;
#endif
        }
        A(i, i) *= temp_coef;
        A(i, i) += penalty + net_charge_constraint;

        // i != j case
        for(int j = i+1; j < n_eff; j++)
        {
            double x_chr_j = pdc.get_x(j);
            double y_chr_j = pdc.get_y(j);
            double z_chr_j = pdc.get_z(j);

            for(int k = 0; k < n_grid_point; k++)
            {
                double x_pot = dx_pot.get_x(x_table[k]);
                double y_pot = dx_pot.get_y(y_table[k]);
                double z_pot = dx_pot.get_z(z_table[k]);
                double disi = std::sqrt((x_pot - x_chr_i)*(x_pot - x_chr_i) +
                                        (y_pot - y_chr_i)*(y_pot - y_chr_i) +
                                        (z_pot - z_chr_i)*(z_pot - z_chr_i));
                double disj = std::sqrt((x_pot - x_chr_j)*(x_pot - x_chr_j) +
                                        (y_pot - y_chr_j)*(y_pot - y_chr_j) +
                                        (z_pot - z_chr_j)*(z_pot - z_chr_j));

#ifdef PDCP_USE_FMATH
                tmp_val = fmath::expd(-invd * (disi + disj)) / (disi * disj);
#else
                tmp_val = std::exp(-invd * (disi + disj)) / (disi * disj);
#endif
                A(i, j) += tmp_val;
            }

            A(i, j) *= temp_coef;
            A(i, j) += net_charge_constraint;
            A(j, i) = A(i, j);
        }
    }

//     std::cout << A << std::endl;

    // Filling vector b
    for(int i = 0; i < n_eff; i++)
    {
        double x_chr_i = pdc.get_x(i);
        double y_chr_i = pdc.get_y(i);
        double z_chr_i = pdc.get_z(i);

        for(int k = 0; k < n_grid_point; k++)
        {
            double x_pot = dx_pot.get_x(x_table[k]);
            double y_pot = dx_pot.get_y(y_table[k]);
            double z_pot = dx_pot.get_z(z_table[k]);
            double disi = std::sqrt((x_pot - x_chr_i)*(x_pot - x_chr_i) +
                                    (y_pot - y_chr_i)*(y_pot - y_chr_i) +
                                    (z_pot - z_chr_i)*(z_pot - z_chr_i));

#ifdef PDCP_USE_FMATH
            tmp_val = fmath::expd(-invd * disi) * phi_table[k] / (disi * diel);
#else
            tmp_val = std::exp(-invd * disi) * phi_table[k] / (disi * diel);
#endif

            b(i) += tmp_val;
        } // k

        b(i) *= factor;
        b(i) += net_charge_constraint * net_charge;
    } // i

    // Solve equation Ax = b
    /////LaLinearSolve(A, x, b);
    Eigen::VectorXd x = A.fullPivLu().solve(b);
    // Output
    std::ofstream ofs(ofname.c_str());
    for(int i = 0; i < n_eff; i++)
    {
        ofs << pdc.get_residue_id(i) << " " << x(i) << std::endl;
    }
    ofs.close();
    return;
}

void calculate_chi(Dx& dx_pot, Pdc& pdc, Charge& charge, Parameter& param)
{
    // Make table
    int i_phi  = 0;
    std::vector<double> phi_table;
    std::vector<int> x_table;
    std::vector<int> y_table;
    std::vector<int> z_table;
    for(int ix = 0; ix < dx_pot.get_grid_x(); ++ix)
    {
        for(int iy = 0; iy < dx_pot.get_grid_y(); ++iy)
        {
            for(int iz = 0; iz < dx_pot.get_grid_z(); ++iz)
            {
                double data = dx_pot.get_data(i_phi);
                if(std::abs(data) > 10e-16)
                {
                    phi_table.push_back(data);
                    x_table.push_back(ix);
                    y_table.push_back(iy);
                    z_table.push_back(iz);
                }
                i_phi++;
            }
        }
    }

    int n_eff = charge.get_n_charge();
    int n_struct = charge.get_n_structure();
    int n_points = phi_table.size();
    // double factor  = 1.87944 * param.get_temp();
    double factor  = 559.219;
    double invd    = 1.0 / param.get_debye();
    double diel    = param.get_diel();

    double tmp_val;

    // Calculate chi
    double chi  = 0.0;
    for(int istruct = 0; istruct < n_struct; ++istruct)
    {
        std::vector<double> phi_deb(n_points, 0.0);

        // Calculate phi_deb
        for(int i = 0; i < n_eff; i++)
        {
            double x_chr_i = pdc.get_x(i);
            double y_chr_i = pdc.get_y(i);
            double z_chr_i = pdc.get_z(i);

            for(int k = 0; k < n_points; k++)
            {
                double x_pot = dx_pot.get_x(x_table[k]);
                double y_pot = dx_pot.get_y(y_table[k]);
                double z_pot = dx_pot.get_z(z_table[k]);
                double disi = std::sqrt((x_pot - x_chr_i)*(x_pot - x_chr_i) +
                                        (y_pot - y_chr_i)*(y_pot - y_chr_i) +
                                        (z_pot - z_chr_i)*(z_pot - z_chr_i));

#ifdef PDCP_USE_FMATH
                tmp_val = factor * charge(istruct, i) * fmath::expd(-invd * disi) / (diel * disi);
#else
                tmp_val = factor * charge(istruct, i) * std::exp(-invd * disi) / (diel * disi);
#endif
                // cout << charge(istruct, i) << "  " << disi << endl;
                phi_deb[k] += tmp_val;
            }
        }

        double chi1 = 0.0;
        double chi2 = 0.0;

        for(int i = 0; i < n_points; ++i)
        {
            double phi = phi_table[i];
            chi1 += (phi_deb[i] - phi)*(phi_deb[i] - phi);
            chi2 += phi * phi;

            if(istruct == 0)
            {
                int tmp_x = x_table[i];
                int tmp_y = y_table[i];
                int grid_center_x = (dx_pot.get_grid_x() - 1) / 2;
                int grid_center_y = (dx_pot.get_grid_y() - 1) / 2;

                if(tmp_x == grid_center_x && tmp_y == grid_center_y)
                {
                    double tmp_z = dx_pot.get_z(z_table[i]);
                    std::cout << tmp_z << " " << phi << " " << phi_deb[i] << std::endl;
                }
            }
        }

        chi += chi1 / chi2;

        if(istruct == 0)
        {
            std::cout << "Chi_native: " << chi << std::endl;
        }
    }

    chi /= static_cast<double>(n_struct);

    std::cout << "Chi: " << chi << std::endl;
}
