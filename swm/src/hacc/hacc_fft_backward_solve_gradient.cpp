/*
 * =====================================================================================
 *
 *       Filename:  hacc_fft_backward_solve_gradient.h
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

#include "hacc_fft_backward_solve_gradient.h"

HaccFFTBackwardSolveGradient::HaccFFTBackwardSolveGradient(
            SWMUserIF* user_if,

            bool* done_from_parent,

            HaccConfig & config,
            double buffer_copy_MBps,
            double fft_work_per_second
            ) :
    HaccFFT(
            user_if,
            done_from_parent,
            config,
            buffer_copy_MBps,
            fft_work_per_second
           )
{}

//BOZO w/ axis
void
HaccFFTBackwardSolveGradient::call() { //backward_solve_gradient(int axis) {
    if (enable_contexts)
    while(1) {
        *done_to_parent = false;
        //kspace_solve_gradient(axis);
        do_fft_compute(2);
        distribution_2_to_3(2);
        distribution_3_to_2(1);
        do_fft_compute(1);
        distribution_2_to_3(1);
        distribution_3_to_2(0);
        do_fft_compute(0);
        distribution_2_to_3(0);
        *done_to_parent = true; yield();
    }
    else
    {
        *done_to_parent = false;
        //kspace_solve_gradient(axis);
        do_fft_compute(2);
        distribution_2_to_3(2);
        distribution_3_to_2(1);
        do_fft_compute(1);
        distribution_2_to_3(1);
        distribution_3_to_2(0);
        do_fft_compute(0);
        distribution_2_to_3(0);
        *done_to_parent = true; yield();

    }
}
