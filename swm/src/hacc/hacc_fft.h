#ifndef _HACC_FFT_H_
#define _HACC_FFT_H_

#include <math.h>

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "app_base_swm_user_code.h"

#include "hacc_config.h"
#include "hacc_ndindex.h"

class SWMProcessUserIF;

typedef enum {
    DISTRIBUTION_3_TO_2,
    DISTRIBUTION_2_TO_3
} direction_t;

class HaccFFT : public SWMUserCode {

    public:

        HaccFFT(
                SWMUserIF* user_if,

                bool* done_from_parent,

                HaccConfig & config,
                double buffer_copy_MBps,
                double fft_work_per_second
                );

        void do_buffer_copy(int chunk_size);

        void do_fft_compute(int /*axis*/);

        void kspace_solve_gradient(int /*axis*/);

        void distribution(int axis, direction_t dir);

        void distribution_3_to_2(int axis);

        void distribution_2_to_3(int axis);


        // Below is code with BLACK MAGIC, directly adapted from HACC

        void pencil_rank_x(const int rank_shape_2d[], const int rank_shape_3d[],
                const int tuple[], int * rank);

        void pencil_rank_y(const int rank_shape_2d[], const int rank_shape_3d[],
                const int tuple[], int * rank);

        void pencil_rank_z(const int rank_shape_2d[], const int rank_shape_3d[],
                const int tuple[], int * rank);

        // More black magic below

        void pencil_tuple_x(const int rank_shape_2d[], const int rank_shape_3d[],
                int rank, int tuple[]);

        void pencil_tuple_y(const int rank_shape_2d[], const int rank_shape_3d[],
                int rank, int tuple[]);

        void pencil_tuple_z(const int rank_shape_2d[], const int rank_shape_3d[],
                int rank, int tuple[]);


    protected:
       
        bool* done_to_parent;
        bool done_to_child;

        HaccConfig & config;
        double buffer_copy_MBps;
        double fft_work_per_second;

        int rank_shape_3d [3];
        int rank_shape_2d [3][3];

        const int LOC2GLOB[3][3] = {{1, 2, 0}, {2, 0, 1}, {0, 1, 2}};

        const int PP_LAXIS[3][2] = {{1, 0}, {1, 0}, {0, 1}};

        const int SIZEOF_ELT = 16; // MPI_DOUBLE_COMPLEX

        RowMajorIndexer<3> indexer_3d;

};

#endif
