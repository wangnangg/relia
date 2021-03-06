#include "helper.h"
#include "easylogging++.h"

template<>
RowSparseM create_matrix<RowSparseM>(uint_t dim, std::vector<double> val)
{
    RowSparseM m(dim);
    uint_t index = 0;
    for (uint_t row_i = 0; row_i < dim; row_i++)
    {
        for (uint_t col_i = 0; col_i < dim; col_i++)
        {
            if (val[index] != 0.0)
            {
                m.add_entry(row_i, col_i, val[index]);
            }
            index++;
        }
    }
    return m;
}

template<>
ColSparseM create_matrix<ColSparseM>(uint_t dim, std::vector<double> val)
{
    ColSparseM m(dim);
    uint_t index = 0;
    for (uint_t row_i = 0; row_i < dim; row_i++)
    {
        for (uint_t col_i = 0; col_i < dim; col_i++)
        {
            if (val[index] != 0.0)
            {
                m.add_entry(col_i, row_i, val[index]);
            }
            index++;
        }
    }
    return m;
}

std::string display(const Vector &v)
{
    std::stringstream ss;
    for (uint_t i = 0; i < v.dim(); i++)
    {
        ss << v(i) << " ";
    }
    ss << std::endl;
    return ss.str();
}
