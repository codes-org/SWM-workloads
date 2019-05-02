
#define __STDC_LIMIT_MACROS
#include <stdint.h>

#include "hacc_timestep.h"

HaccTimestep::HaccTimestep (
        SWMUserIF* user_if,

        bool* done_from_parent,

        HaccConfig & config,

        const double ninteractions_per_rank_mean,
        const double ninteractions_per_rank_delta,
        const double ninteractions_per_rank_per_wallsecond,

        const double buffer_copy_MBps,
        const double fft_work_per_second,

        const bool enable_hacc_fft,
        const bool enable_hacc_exchange,
        const bool enable_hacc_checksum

        /*
        HaccComputeRCBTree & rcb,
        HaccFFTBackwardSolveGradient & fft_backward,
        HaccFFTForwardSolve & fft_forward,
        HaccExchange & exchange
        */
        ) :
    SWMUserCode(user_if),
    done_to_parent(done_from_parent),
    config(config),
    enable_hacc_fft(enable_hacc_fft),
    enable_hacc_exchange(enable_hacc_exchange),
    enable_hacc_checksum(enable_hacc_checksum)
{ 

    // Compute models
    rcb = new HaccComputeRCBTree (
            user_if,
            &done_to_child,
            config,
            ninteractions_per_rank_mean,
            ninteractions_per_rank_delta,
            ninteractions_per_rank_per_wallsecond
            );

    fft_forward_solve = new HaccFFTForwardSolve (
            user_if,
            &done_to_child,
            config,
            buffer_copy_MBps,
            fft_work_per_second
            );

    fft_backward_solve_gradient = new HaccFFTBackwardSolveGradient (
            user_if,
            &done_to_child,
            config,
            buffer_copy_MBps,
            fft_work_per_second
            );

    exchange = new HaccExchange (
            user_if,
            &done_to_child,
            config,
            buffer_copy_MBps
            );


    /*
    cerr << "Hacc sizes: " << endl;
    cerr << "Hacc pid: " << process->process_id << endl;
    cerr << "Hacc timestep: " << sizeof(HaccTimestep) << endl;
    cerr << "Hacc rcb: " << sizeof(HaccComputeRCBTree) << endl;
    cerr << "Hacc fft_forward_solve: " << sizeof(HaccFFTForwardSolve) << endl;
    cerr << "Hacc fft_backward_solve: " << sizeof(HaccFFTBackwardSolveGradient) << endl;
    cerr << "Hacc exchange : " << sizeof(HaccExchange) << endl;
    */
}

void
HaccTimestep::sub_cycle() {
    if (enable_contexts)
    while(1) {
        (*rcb)(); //rcb.build_tree_and_evaluate_forces();
        if(done_to_child) break;
        else yield();
    }
    else
        (*rcb)(); //rcb.build_tree_and_evaluate_forces();
}

void
HaccTimestep::map2_poisson_forward() {
    // TODO:compute: forward CIC
    if (enable_hacc_fft)
    {
        if (enable_contexts)
        while(1) {
            (*fft_forward_solve)(); //fft.forward_solve();
            if(done_to_child) break;
            else yield();
        }
        else
            (*fft_forward_solve)(); //fft.forward_solve();
    }
}

void
HaccTimestep::map2_poisson_backward_gradient() {
    for (int idim=0; idim<3; idim++) {
        if (enable_hacc_fft)
        {
            if (enable_contexts)
            while(1) {
                (*fft_backward_solve_gradient)(); //fft.backward_solve_gradient(idim);
                if(done_to_child) break;
                else yield();
            }
            else
                (*fft_backward_solve_gradient)(); //fft.backward_solve_gradient(idim);
        }

        if (enable_hacc_exchange)
        {
            // TODO:compute: FFT copy array backwards
            if (enable_contexts)
            while(1) {
                (*exchange)();
                if(done_to_child) break;
                else yield();
            }
            else
                (*exchange)();
            // TODO:compute: inverse CIC
        }
    }
}


void
HaccTimestep::call() { //do_steps() {
    *done_to_parent = false;

#if 0
    //use this just to test nested contexts...
    for (int istep=0; istep<2; istep++) {
        SWM_Compute(2);
        std::cout << "random is " << user_if->RNGIF()->Get(100) << std::endl;
    }
#endif

//#if 0
    for (int istep=0; istep<nstep; istep++) {

        // Half-kick at first step
        if (istep == 0) {
            map2_poisson_forward();
            map2_poisson_backward_gradient();
        }

        // Sub-stepping
        for (int isubstep=0;isubstep<nsub; isubstep++) {
            sub_cycle();
        }

        // FFT memory reallocation
        if (enable_hacc_fft && do_drop_memory) {
            // FFT memory reallocation will rebuild the solver, which calls
            // MPI_Cart_create.  MPI_Cart_create is a collective op which
            // contains an implicit synchronization and acts as a barrier,
            // so it's important to emulate its effect
            //backend.comm_emulate_cart_create();
            SWM_Barrier(
                    SWM_COMM_WORLD,  //comm_id
                    config.request_vc,  //reqvc
                    config.response_vc  //rspvc
                    );
        }

        // Checksum particles
        //backend.comm_allreduce(8); // reduce a single MPI_LONG_LONG
        if (enable_hacc_checksum)
        SWM_Allreduce(
                8,  //msg_size
                config.pkt_rsp_bytes, // pkt_rsp_bytes
                SWM_COMM_WORLD,  //comm_id
                config.request_vc,  //reqvc
                config.response_vc,  //rspvc
                NO_BUFFER,  //sendbuf
                NO_BUFFER   //recvbuf
                );
        
        // Get rho into spectral domain
        map2_poisson_forward();

        // Checksum grid density
        //backend.comm_allreduce(8); // reduce a single MPI_DOUBLE
        if (enable_hacc_checksum)
         SWM_Allreduce(
                8,  //msg_size
                config.pkt_rsp_bytes, // pkt_rsp_bytes
                SWM_COMM_WORLD,  //comm_id
                config.request_vc,  //reqvc
                config.response_vc,  //rspvc
                NO_BUFFER,  //sendbuf
                NO_BUFFER   //recvbuf
                );
       
        // P(k) computation may occur here, and does a couple of Allreduce
        // on arrays of doubles It seems P(k) is never computed within the
        // timesteps though, and the .ini and .fin power spectra are
        // computed outside of the time loop.
        
        // Solve gradient and back to real domain
        map2_poisson_backward_gradient();

        // AFAIK, this is never called
        // map2_poisson_backward_potential();

        if (enable_hacc_checksum)
        {
            for(int ar=0; ar<3; ar++) {
                SWM_Allreduce(
                        8,  //msg_size
                        config.pkt_rsp_bytes, // pkt_rsp_bytes
                        SWM_COMM_WORLD,  //comm_id
                        config.request_vc,  //reqvc
                        config.response_vc,  //rspvc
                        NO_BUFFER,  //sendbuf
                        NO_BUFFER   //recvbuf
                        );
            }
            SWM_Barrier(
                    SWM_COMM_WORLD,  //comm_id
                    config.request_vc,  //reqvc
                    config.response_vc  //rspvc
                    );
        }
    }

    *done_to_parent = true; yield();
}
