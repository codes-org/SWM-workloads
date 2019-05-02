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

#ifndef _HACC_FFT_BACKWARD_SOLVE_GRADIENT_H_
#define _HACC_FFT_BACKWARD_SOLVE_GRADIENT_H_

#include "hacc_fft.h"

class HaccFFT;

class HaccFFTBackwardSolveGradient : public HaccFFT {

    public:

        HaccFFTBackwardSolveGradient(
                    SWMUserIF* user_if,

                    bool* done_from_parent,

                    HaccConfig & config,
                    double buffer_copy_MBps,
                    double fft_work_per_second
                    );

        //BOZO w/ axis
        void call();

};

#endif
