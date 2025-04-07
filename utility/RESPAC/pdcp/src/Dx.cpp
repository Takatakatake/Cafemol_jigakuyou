#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <boost/algorithm/string.hpp>
#include "Dx.hpp"

Dx::Dx(const std::string& fname)
{
    std::ifstream ifs(fname.c_str());
    if(!ifs.good())
    {
        std::cerr << "The file \"" << fname << "\" cannot open" << std::endl;
        throw std::invalid_argument("file open error");
    }
 
    std::string line;
    std::vector<std::string> columns;
    std::size_t idata(0);
    int num_mdelta_read(1);
 
    while(std::getline(ifs, line))
    {
        if(line.substr(0, 8) == "object 1")
        {
            boost::algorithm::split(
                        columns,
                        line,
                        boost::is_any_of(" "),
                        boost::token_compress_on
                    );
            m_grid.at(0) = std::stoi(columns[5]);
            m_grid.at(1) = std::stoi(columns[6]);
            m_grid.at(2) = std::stoi(columns[7]);
        }
        else if(line.substr(0, 6) == "origin")
        {
            boost::algorithm::split(
                        columns,
                        line,
                        boost::is_any_of(" "),
                        boost::token_compress_on
                    );
            m_origin.at(0) = std::stod(columns[1]);
            m_origin.at(1) = std::stod(columns[2]);
            m_origin.at(2) = std::stod(columns[3]);
        }
        else if(line.substr(0, 5) == "delta")
        {
            boost::algorithm::split(
                        columns,
                        line,
                        boost::is_any_of(" "),
                        boost::token_compress_on
                    );
            // XXX: read diagonal value only
            m_delta.at(num_mdelta_read - 1)
                = std::stod(columns[num_mdelta_read]);
            ++num_mdelta_read;
        }
        else if(line.substr(0, 8) == "object 3")
        {
            boost::algorithm::split(
                        columns,
                        line,
                        boost::is_any_of(" "),
                        boost::token_compress_on
                    );
            m_n_data = std::stoi(columns[9]);
            m_data.resize(m_n_data);
        }
        else if(line.find_first_of("-0123456789") == 0)
        {
            boost::algorithm::split(
                        columns,
                        line,
                        boost::is_any_of(" "),
                        boost::algorithm::token_compress_on
                    );
            for(auto iter = columns.begin(); iter != columns.end(); ++iter)
            {
                // each line has a space at end of the line
                if(iter->empty()) continue;
                m_data.at(idata) = std::stod(*iter);
                idata++;
                /* XXX: when the number of data(m_n_data) is not a multiple  *
                 *      of 3, the original code writes a empty value in      *
                 *      the out-of-range element. In the case, at(size_type) *
                 *      method throws out-of-range exception.                *
                 *      the reason why not m_data.at(idata++) is readability */
            }
        }
    }
    assert(idata == m_n_data);
 
    for(std::size_t i = 0; i < 3; i++)
    {
        m_grid_center.at(i) =
            m_origin.at(i) + 0.5 * (m_grid.at(i) - 1) * m_delta.at(i);
    }
 
    for(std::size_t dim_index = 0; dim_index < 3; ++dim_index)
    {
        const double offset =
            -0.5 * (m_grid.at(dim_index) - 1) * m_delta.at(dim_index) +
             m_grid_center.at(dim_index);

        const double delta_length = m_delta.at(dim_index);

        for(std::size_t grid_index = 0;
            grid_index < m_grid.at(dim_index);
            ++grid_index)
        {
            m_coord.at(dim_index).push_back(
                    offset + grid_index * delta_length
                );
        }
    }
//                CPU                   :
    std::cerr << "Info >> Dx" << std::endl;
    std::cerr << "The # of Grid         : " << m_grid[0]        << " "
                                            << m_grid[1]        << " " 
                                            << m_grid[2]        << std::endl;
    std::cerr << "Origin                : " << m_origin[0]      << " "
                                            << m_origin[1]      << " "
                                            << m_origin[2]      << std::endl;
    std::cerr << "Grid Center           : " << m_grid_center[0] << " "
                                            << m_grid_center[1] << " "
                                            << m_grid_center[2] << std::endl;
    std::cerr << "Delta                 : " << m_delta[0]       << " "
                                            << m_delta[1]       << " " 
                                            << m_delta[2]       << std::endl;
    std::cerr << "The # of Data         : " << m_n_data         << std::endl;
    std::cerr << std::endl;

    return;
}

Dx Dx::operator*(const Dx& other)
{
    Dx result = *this;
    // int b = 0;
    for(std::size_t i = 0; i < other.get_n_data(); ++i)
    {
        // if (result.get_data(i) <=0.1 && other.get_data(i) >=0.9 && b++ < 100) {
        //     cout << i << "  " << result.get_data(i) << endl;
        // }
        result.set_data(i, result.get_data(i) * other.get_data(i));
    }
    return result;
}
