#ifndef _HACC_EXCHANGE_H_
#define _HACC_EXCHANGE_H_

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "app_base_swm_user_code.h"

#include "hacc_config.h"
#include "hacc_ndindex.h"

class SWMProcessUserIF;

class HaccExchange : public SWMUserCode {

    protected:

    bool* done_to_parent;

    HaccConfig & config;
    double buffer_copy_MBps;

    RowMajorIndexer<3> index3d;

    int buffer_size;
    int ng_alive[3];
    int mytuple[3];

    const int SIZEOF_ELT = 4; // We have -DGRID_32, so we're using MPI_FLOAT
    static const int NUM_OF_NEIGHBORS = 26;

    int neighbor_rank[NUM_OF_NEIGHBORS];

    public:

        HaccExchange(
                SWMUserIF* user_if,

                bool* done_from_parent,

                HaccConfig & config,
                double buffer_copy_MBps
                );

        int get_my_neighbor(int ishift, int jshift, int kshift);

        void do_buffer_copy(int buffer_size);

        void exchange(int inbor_send_to, int inbor_recv_from);

        void call();

};


#endif
