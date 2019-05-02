/*
 * =====================================================================================
 *
 *       Filename:  hacc_fft_forward_solve.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  08/06/2013 12:53:36 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  John Thompson (FPD), john.d.thompson@intel.com
 *        Company:  Intel
 *
 * =====================================================================================
 */

#include "hacc_fft_forward_solve.h"

HaccFFTForwardSolve::HaccFFTForwardSolve(
            SWMUserIF* user_if,

            bool* done,

            HaccConfig & config,
            double buffer_copy_MBps,
            double fft_work_per_second
            ) :
    HaccFFT(
            user_if,
            done,
            config,
            buffer_copy_MBps,
            fft_work_per_second
           )
{}

void
HaccFFTForwardSolve::call() { //forward_solve() {

    if (enable_contexts)
    while(1) {
        *done_to_parent = false;
        distribution_3_to_2(0);
        do_fft_compute(0);
        distribution_2_to_3(0);
        distribution_3_to_2(1);
        do_fft_compute(1);
        distribution_2_to_3(1);
        distribution_3_to_2(2);
        do_fft_compute(2);
        *done_to_parent = true; yield();
    }
    else
    {
        *done_to_parent = false;
        distribution_3_to_2(0);
        do_fft_compute(0);
        distribution_2_to_3(0);
        distribution_3_to_2(1);
        do_fft_compute(1);
        distribution_2_to_3(1);
        distribution_3_to_2(2);
        do_fft_compute(2);
        *done_to_parent = true; yield();

    }
}
