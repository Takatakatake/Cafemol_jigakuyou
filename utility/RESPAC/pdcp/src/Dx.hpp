#ifndef RESPAC_PDCP_DX_H
#define RESPAC_PDCP_DX_H

#include <string>
#include <vector>
#include <array>

class Dx
{
    public:

        Dx(const std::string& fname);
        Dx operator*(const Dx& other);

        double get_data(const int id) const {return m_data.at(id);}
        void   set_data(const int id, const double data){m_data.at(id) = data;}

        std::size_t get_n_data() const {return m_n_data;}

        int get_grid_x() const {return m_grid[0];}
        int get_grid_y() const {return m_grid[1];}
        int get_grid_z() const {return m_grid[2];}

        double get_x(const int id) const {return m_coord.at(0).at(id);}
        double get_y(const int id) const {return m_coord.at(1).at(id);}
        double get_z(const int id) const {return m_coord.at(2).at(id);}

    private:

        // m_data size
        std::size_t m_n_data;
        std::array<int, 3> m_grid;
        std::array<double, 3> m_grid_center;
        std::array<double, 3> m_origin;
        std::array<double, 3> m_delta;

        // .at(dimention_index).at(grid_index)
        std::array<std::vector<double>, 3> m_coord;

        // core data
        // size of this vector is same as m_n_data;
        std::vector<double> m_data;

};

#endif //RESPAC_PDCP_DX_H
