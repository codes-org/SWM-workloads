
#include "hacc_config.h"

HaccConfig::HaccConfig (
        int _ng, double _box_length,
        int _nranks, int _myrank,
        const int * _rank_shape_3d,
        const int * _rank_shape_2d_x,
        const int * _rank_shape_2d_y,
        const int * _rank_shape_2d_z,
        uint32_t request_vc,
        uint32_t response_vc,
        uint32_t pkt_rsp_bytes) :
    ng(_ng),
    box_length(_box_length),
    nranks(_nranks),
    myrank(_myrank),
    request_vc(request_vc),
    response_vc(response_vc),
    pkt_rsp_bytes(pkt_rsp_bytes)
{
    for (int i=0; i<NDIM; i++) {
        rank_shape_3d  [i] = _rank_shape_3d  [i];
        rank_shape_2d_x[i] = _rank_shape_2d_x[i];
        rank_shape_2d_y[i] = _rank_shape_2d_y[i];
        rank_shape_2d_z[i] = _rank_shape_2d_z[i];
    }

    rank_shape_1d[0] = _myrank;
    rank_shape_1d[1] = 1;
    rank_shape_1d[2] = 1;

    overload_len = 8.0; // Fixed in CORAL indat file
    phys2grid = ng / box_length;
    ng_overload = static_cast< int > (ceilf(overload_len*phys2grid));
}
