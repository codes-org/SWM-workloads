
#include "hacc_fft.h"
//#include <asim/restricted/swm_process.h>

HaccFFT::HaccFFT(
        SWMUserIF* user_if,

        bool* done_from_parent,

        HaccConfig & config,
        double buffer_copy_MBps,
        double fft_work_per_second
        ) :
    SWMUserCode(user_if),
    done_to_parent(done_from_parent),
    config(config),
    buffer_copy_MBps(buffer_copy_MBps),
    fft_work_per_second(fft_work_per_second)//,
    //indexer_3d()
{

    //BOZO -- this should be done in some common config?
    for (int i=0; i<3; i++) {
        rank_shape_3d[i] = config.rank_shape_3d[i];

        rank_shape_2d[0][i] = config.rank_shape_2d_x[i];
        rank_shape_2d[1][i] = config.rank_shape_2d_y[i];
        rank_shape_2d[2][i] = config.rank_shape_2d_z[i];
    }

    // Sanity checks
    for (int i=0; i<3; i++) {
        // Check 3D decomp divisibility conditions
        assert(config.ng % rank_shape_3d[i] == 0);

        // Check 2D decomp divisibility conditions
        assert(config.ng % rank_shape_2d[0][i] == 0);
        assert(config.ng % rank_shape_2d[1][i] == 0);
        assert(config.ng % rank_shape_2d[2][i] == 0);
    }
    for (int axis=0; axis<3; axis++) {
        // Check pencil decomp embedding conditions
        const int uax = LOC2GLOB[axis][0];
        const int vax = LOC2GLOB[axis][1];

        const int ng2u = config.ng / rank_shape_2d[axis][uax];
        const int ng2v = config.ng / rank_shape_2d[axis][vax];

        const int ng3u = config.ng / rank_shape_3d[uax];
        const int ng3v = config.ng / rank_shape_3d[vax];

        assert(ng3u % ng2u == 0);
        assert(ng3v % ng2v == 0);
    }
}

void
HaccFFT::do_buffer_copy(int chunk_size) 
{
    const double seconds = (chunk_size*SIZEOF_ELT/1e6) / buffer_copy_MBps;
    SWM_Compute(seconds);
}

void
HaccFFT::do_fft_compute(int /*axis*/) {
    // Number of FFTs for one rank
    const long nfft_per_rank = ((long)config.ng)*((long)config.ng) / config.nranks;

    // Length of one FFT
    const long fft_size = config.ng;

    // Work associated with the FFTs for one rank
    const double fft_work = nfft_per_rank * fft_size*log(fft_size);

    const double seconds = fft_work / fft_work_per_second;
    SWM_Compute(seconds);
}

void
HaccFFT::kspace_solve_gradient(int /*axis*/) {
    // TODO:compute: kspace_solve_gradient
}

void
HaccFFT::distribution(int axis, direction_t dir) {
    assert(axis >=0 && axis < 3);

    int block_shape_2d[3];
    int block_shape_3d[3];
    int chunk_shape[3];

    for (int i=0; i<3; i++) {
        block_shape_2d[i] = config.ng / rank_shape_2d[axis][i];
        block_shape_3d[i] = config.ng / rank_shape_3d[i];

        chunk_shape[i] = block_shape_2d[i];
    }
    chunk_shape[axis] = config.ng / rank_shape_3d[axis];

    long chunk_size = 1;
    for (int i=0; i<3; i++) chunk_size *= chunk_shape[i];

    // Coordinates (tuples) of this rank in 2D and 3D decomps
    int self_tup_2d[3];
    int self_tup_3d[3];
    switch (axis) {
        case 0:
            pencil_tuple_x(rank_shape_2d[axis], rank_shape_3d, config.myrank, self_tup_2d); break;
        case 1:
            pencil_tuple_y(rank_shape_2d[axis], rank_shape_3d, config.myrank, self_tup_2d); break;
        case 2:
            pencil_tuple_z(rank_shape_2d[axis], rank_shape_3d, config.myrank, self_tup_2d); break;
        default:
            abort(); break;
    }
    indexer_3d.index_to_tuple(rank_shape_3d, config.myrank, self_tup_3d);
    
    const int npeers = config.rank_shape_3d[axis];

    // Local counter for chunks within a 3D cube
    int pp[2] = {0, 0};
    const int p1_gaxis = LOC2GLOB[axis][PP_LAXIS[axis][1]];
    const int p1max = rank_shape_2d[axis][p1_gaxis]/rank_shape_3d[p1_gaxis] - 1;

    for (int p=0; p<npeers; p++) {

        int grid_coords[3];

        // Compute the tuple of the current chunk (peer) in the 3D decomp
        // Consecutive p's are traversing the pencil length along `axis`
        for (int i=0; i<3; i++) {
            grid_coords[i] = self_tup_2d[i] * block_shape_2d[i];
        }
        grid_coords[axis] += p * block_shape_3d[axis];

        int p_tup_3d[3];
        for (int i=0; i<3; i++) {
            p_tup_3d[i] = grid_coords[i] / block_shape_3d[i];
        }

        // Compute the tuple of the current chunk (peer) in the 2D decomp
        for (int i=0; i<3; i++) {
            grid_coords[i] = self_tup_3d[i] * block_shape_3d[i];
        }
        for (int i=0; i<2; i++) { // Loop bound is really TWO here
            const int iglob = LOC2GLOB[axis][i];
            grid_coords[iglob] += pp[PP_LAXIS[axis][i]] * block_shape_2d[iglob];
        }
        // Increment pp counter
        if (pp[1] == p1max) { pp[0]++; pp[1] = 0; } else { pp[1]++; }

        int p_tup_2d[3];
        for (int i=0; i<3; i++) {
            p_tup_2d[i] = grid_coords[i] / block_shape_2d[i];
        }
        p_tup_2d[axis] = 0;

        // Now compute the corresponding rank ids
        int p_rank_2d;
        int p_rank_3d;
        switch (axis) {
            case 0:
                pencil_rank_x(rank_shape_2d[axis], rank_shape_3d, p_tup_2d, &p_rank_2d); break;
            case 1:
                pencil_rank_y(rank_shape_2d[axis], rank_shape_3d, p_tup_2d, &p_rank_2d); break;
            case 2:
                pencil_rank_z(rank_shape_2d[axis], rank_shape_3d, p_tup_2d, &p_rank_2d); break;
            default:
                abort(); break;
        }
        indexer_3d.tuple_to_index(rank_shape_3d, p_tup_3d, &p_rank_3d);

        int send_peer, recv_peer;
        switch (dir) {
            case DISTRIBUTION_2_TO_3:
                recv_peer = p_rank_3d;
                send_peer = p_rank_2d;
                break;
            case DISTRIBUTION_3_TO_2:
                recv_peer = p_rank_2d;
                send_peer = p_rank_3d;
                break;
            default:
                abort();
                break;
        }

        // Pack data into send buffers
        // HACC has different cases here based on distribution direction.
        // Ignoring the potential impact of buffer shape on perf, all cases
        // will perform a copy of size chunk_size. In addition, each rank will
        // perform exactly one such copy here.
        do_buffer_copy(chunk_size);

        // Perform the MPI exchange
        // FFT part of comm parttern
        SWM_Sendrecv(
                SWM_COMM_WORLD, //0,  //comm_id
                send_peer,  //sendpeer
                0,  //sendtag
                config.request_vc,  //sendreqvc
                config.response_vc,  //sendrspvc
                NO_BUFFER,  //sendbuf
                chunk_size*SIZEOF_ELT,  //sendbytes
                config.pkt_rsp_bytes,   //rspbytes
                recv_peer,  //recvpeer
                0,  //recvtag
                NO_BUFFER   //recvbuf
                );


        // Unpack data from receive buffers
        // See note in comment above
        do_buffer_copy(chunk_size);

    }

}

void
HaccFFT::distribution_3_to_2(int axis) {
    distribution(axis, DISTRIBUTION_3_TO_2);
}

void
HaccFFT::distribution_2_to_3(int axis) {
    distribution(axis, DISTRIBUTION_2_TO_3);
}


// Below is code with BLACK MAGIC, directly adapted from HACC

void
HaccFFT::pencil_rank_x(const int rank_shape_2d[], const int rank_shape_3d[],
        const int tuple[], int * rank) {

    int num_pen_in_cube_col = rank_shape_2d[1] / rank_shape_3d[1];
    int num_pen_in_cube_row = rank_shape_2d[2] / rank_shape_3d[2];
    assert(num_pen_in_cube_col != 0 && num_pen_in_cube_row != 0);

    int alpha = tuple[1] % num_pen_in_cube_col;
    int beta  = tuple[2] % num_pen_in_cube_row;
    int num_cubes = rank_shape_3d[2] * rank_shape_3d[1];

    *rank = (alpha*num_cubes)
        + ((tuple[1]/num_pen_in_cube_col)*rank_shape_3d[2])
        + (beta*(num_cubes)*num_pen_in_cube_col) + tuple[2]/num_pen_in_cube_row;
}

void
HaccFFT::pencil_rank_y(const int rank_shape_2d[], const int rank_shape_3d[],
        const int tuple[], int * rank) {

    int num_pen_in_cube_col = rank_shape_2d[0] / rank_shape_3d[0];
    int num_pen_in_cube_row = rank_shape_2d[2] / rank_shape_3d[2];
    assert(num_pen_in_cube_col != 0 && num_pen_in_cube_row != 0);

    int beta = tuple[2] % num_pen_in_cube_row;
    *rank = tuple[0] * rank_shape_2d[2]
        + beta * rank_shape_3d[2]
        + tuple[2] / num_pen_in_cube_row;
}

void
HaccFFT::pencil_rank_z(const int rank_shape_2d[], const int rank_shape_3d[],
        const int tuple[], int * rank) {

    int num_pen_in_cube_col = rank_shape_2d[1] / rank_shape_3d[1];
    int num_pen_in_cube_row = rank_shape_2d[0] / rank_shape_3d[0];
    int num_pen_in_cube = rank_shape_3d[2];
    assert(num_pen_in_cube_col != 0 && num_pen_in_cube_row != 0);

    int alpha = tuple[1] % num_pen_in_cube_col;
    int beta  = tuple[0] % num_pen_in_cube_row;
    *rank = alpha
        + ((tuple[1]/num_pen_in_cube_col)*num_pen_in_cube)
        + (beta*num_pen_in_cube_col)
        + (tuple[0]/num_pen_in_cube_row)*rank_shape_2d[1]*num_pen_in_cube_row;
}

// More black magic below

void
HaccFFT::pencil_tuple_x(const int rank_shape_2d[], const int rank_shape_3d[],
        int rank, int tuple[]) {

    assert(rank_shape_2d[0] == 1);
    tuple[0] = 0;
    int num_pen_in_cube_col = rank_shape_2d[1] / rank_shape_3d[1];
    int num_pen_in_cube_row = rank_shape_2d[2] / rank_shape_3d[2];
    int num_cubes = (rank_shape_3d[2]*rank_shape_3d[1]);

    int num_repeats=rank/(num_cubes);
    int low_rank=rank-num_repeats*num_cubes;
    tuple[1] = (low_rank/rank_shape_3d[2])*num_pen_in_cube_col
        + num_repeats%num_pen_in_cube_col;
    tuple[2] = (low_rank%rank_shape_3d[2])*num_pen_in_cube_row + num_repeats/num_pen_in_cube_col;

}

void
HaccFFT::pencil_tuple_y(const int rank_shape_2d[], const int rank_shape_3d[],
        int rank, int tuple[]) {

    assert(rank_shape_2d[1] == 1);
    tuple[1] = 0;
    int num_pen_in_cube_row = rank_shape_2d[2]/rank_shape_3d[2];
    int alpha = rank%(rank_shape_2d[2]);
    tuple[0] = rank/rank_shape_2d[2];

    tuple[2] = (alpha/rank_shape_3d[2])
        + (alpha%rank_shape_3d[2])*num_pen_in_cube_row;
}

void
HaccFFT::pencil_tuple_z(const int rank_shape_2d[], const int rank_shape_3d[],
        int rank, int tuple[]) {

    assert(rank_shape_2d[2] == 1);
    tuple[2] = 0;
    int num_pen_in_cube_col = rank_shape_2d[1]/rank_shape_3d[1];
    int num_pen_in_cube_row = rank_shape_2d[0]/rank_shape_3d[0];
    int num_pen_in_cube = rank_shape_3d[2];
    int alpha = rank/(rank_shape_2d[1]*num_pen_in_cube_row);
    tuple[0] = alpha*num_pen_in_cube_row + (rank%num_pen_in_cube)/num_pen_in_cube_col;
    tuple[1] =
        ((rank%(rank_shape_2d[1]*num_pen_in_cube_row))/num_pen_in_cube)*num_pen_in_cube_col + rank%num_pen_in_cube_col;
}
