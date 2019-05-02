
#include "hacc_exchange.h"

HaccExchange::HaccExchange(
        SWMUserIF* user_if,

        bool* done_from_parent,

        HaccConfig & config,
        double buffer_copy_MBps
        ) :
    SWMUserCode(user_if),
    done_to_parent(done_from_parent),
    config(config),
    buffer_copy_MBps(buffer_copy_MBps),
    index3d()
{
    
    // Set our position in the 3D rank layout
    index3d.index_to_tuple(config.rank_shape_3d, config.myrank, mytuple);

    //const int dead0 = config.ng_overload;
    //const int dead1 = config.ng_overload+1;

    // Compute buffer size (max of all halo sizes)
    int max_ng = 0;
    for (int i=0; i<3; i++) {
        ng_alive[i] = config.ng/config.rank_shape_3d[i];
        if (ng_alive[i] > max_ng) max_ng = ng_alive[i];
    }

    //Thomas needs to check this change for a compiler warning... BOZO, JOHNT
    //buffer_size = max_ng * max_ng * (dead0 > dead1 ? dead0 : dead1);
    buffer_size = max_ng * max_ng * (config.ng_overload == INT_MAX ? config.ng_overload : (config.ng_overload+1));

    // Fill in neighbor list
    // CAUTION: Do not change the order of the following lines
    // unless you know what you're doing!
    int i = 0;

    // Face neighbors
    neighbor_rank[i] = get_my_neighbor(-1,  0,  0); i++;
    neighbor_rank[i] = get_my_neighbor( 1,  0,  0); i++;
    neighbor_rank[i] = get_my_neighbor( 0, -1,  0); i++;
    neighbor_rank[i] = get_my_neighbor( 0,  1,  0); i++;
    neighbor_rank[i] = get_my_neighbor( 0,  0, -1); i++;
    neighbor_rank[i] = get_my_neighbor( 0,  0,  1); i++;

    // Edge neighbors
    neighbor_rank[i] = get_my_neighbor(-1, -1,  0); i++;
    neighbor_rank[i] = get_my_neighbor( 1,  1,  0); i++;
    neighbor_rank[i] = get_my_neighbor(-1,  1,  0); i++;
    neighbor_rank[i] = get_my_neighbor( 1, -1,  0); i++;

    neighbor_rank[i] = get_my_neighbor( 0, -1, -1); i++;
    neighbor_rank[i] = get_my_neighbor( 0,  1,  1); i++;
    neighbor_rank[i] = get_my_neighbor( 0, -1,  1); i++;
    neighbor_rank[i] = get_my_neighbor( 0,  1, -1); i++;

    neighbor_rank[i] = get_my_neighbor(-1,  0, -1); i++;
    neighbor_rank[i] = get_my_neighbor( 1,  0,  1); i++;
    neighbor_rank[i] = get_my_neighbor( 1,  0, -1); i++;
    neighbor_rank[i] = get_my_neighbor(-1,  0,  1); i++;

    // Corner neighbors
    neighbor_rank[i] = get_my_neighbor(-1, -1, -1); i++;
    neighbor_rank[i] = get_my_neighbor( 1,  1,  1); i++;
    neighbor_rank[i] = get_my_neighbor(-1, -1,  1); i++;
    neighbor_rank[i] = get_my_neighbor( 1,  1, -1); i++;
    neighbor_rank[i] = get_my_neighbor(-1,  1, -1); i++;
    neighbor_rank[i] = get_my_neighbor( 1, -1,  1); i++;
    neighbor_rank[i] = get_my_neighbor(-1,  1,  1); i++;
    neighbor_rank[i] = get_my_neighbor( 1, -1, -1); i++;

    assert(i == NUM_OF_NEIGHBORS);
}

int
HaccExchange::get_my_neighbor(int ishift, int jshift, int kshift) {
    int nbtuple[3] = { mytuple[0]+ishift, mytuple[1]+jshift, mytuple[2]+kshift };
    // Handle periodicity
    for (int i=0; i<3; i++) {
        nbtuple[i] = (nbtuple[i]+config.rank_shape_3d[i]) % config.rank_shape_3d[i];
    }
    int nbrank;
    index3d.tuple_to_index(config.rank_shape_3d, nbtuple, &nbrank);
    return nbrank;
}

void
HaccExchange::do_buffer_copy(int buffer_size) {
    const int nbytes = buffer_size * SIZEOF_ELT;
    const double seconds = (nbytes/1e6) / buffer_copy_MBps;
    SWM_Compute(seconds);
}

void
HaccExchange::exchange(int inbor_send_to, int inbor_recv_from) {
    int rank_send_to   = neighbor_rank[inbor_send_to  ];
    int rank_recv_from = neighbor_rank[inbor_recv_from];
    int nbytes = buffer_size * SIZEOF_ELT;
    uint32_t rsp_bytes = 0; // should really be nonzero to have network send rsp pkt when matching completes

    // Important note about buffer copies:
    // While the MPI op size is exactly buffer_size, the copy size is not (it
    // can be smaller).  However, because the buffer copies are not a major
    // hotspot, we overlook this fact for simplicity and use buffer_size for
    // the copies as well. The idea is just to introduce some delay between the
    // MPI comms.

    // Pack send buffer
    do_buffer_copy(buffer_size);

    // Perform the MPI exchange
    SWM_Sendrecv(
            SWM_COMM_WORLD, //0,  //comm_id
            rank_send_to,  //sendpeer
            0,  //sendtag
            config.request_vc,  //sendreqvc
            config.response_vc,  //sendrspvc
            NO_BUFFER,  //sendbuf
            nbytes,  //sendbytes
            rsp_bytes, //rspbytes
            rank_recv_from,  //recvpeer
            0,  //recvtag
            NO_BUFFER   //recvbuf
            );



    // Unpack receive buffer
    do_buffer_copy(buffer_size);
}

void
HaccExchange::call() { //exchange_grid() {
    if(enable_contexts)
    while(1) {
        *done_to_parent = false;
        for (int inbor=0; inbor<NUM_OF_NEIGHBORS; inbor+=2) {
            exchange(inbor, inbor+1);
            exchange(inbor+1, inbor);
        }
        *done_to_parent = true; yield();
    }
    else
    {
        *done_to_parent = false;
        for (int inbor=0; inbor<NUM_OF_NEIGHBORS; inbor+=2) {
            exchange(inbor, inbor+1);
            exchange(inbor+1, inbor);
        }
        *done_to_parent = true; yield();
    }
}
