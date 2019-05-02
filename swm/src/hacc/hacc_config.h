#ifndef _HACC_CONFIG_HPP
#define _HACC_CONFIG_HPP

#include <stdint.h>
#include <math.h>

const int NDIM = 3;

class HaccConfig {
    public:

        int ng;
        double box_length;

        int nranks;
        int myrank;

        int rank_shape_3d  [NDIM];

        int rank_shape_2d_x[NDIM];
        int rank_shape_2d_y[NDIM];
        int rank_shape_2d_z[NDIM];

        int rank_shape_1d  [NDIM];

        double overload_len;
        double phys2grid;

        int ng_overload;

        uint32_t request_vc;
        uint32_t response_vc;
        uint32_t pkt_rsp_bytes;

        HaccConfig (
                int _ng, double _box_length,
                int _nranks, int _myrank,
                const int * _rank_shape_3d,
                const int * _rank_shape_2d_x,
                const int * _rank_shape_2d_y,
                const int * _rank_shape_2d_z,
                uint32_t request_vc,
                uint32_t response_vc,
                uint32_t pkt_rsp_bytes);

};


#endif
