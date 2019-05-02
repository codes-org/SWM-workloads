/*
 * =====================================================================================
 *
 *       Filename:  hacc_fft_forware_solve.h
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

#ifndef _HACC_FFT_FORWARD_SOLVE_H_
#define _HACC_FFT_FORWARD_SOLVE_H_

#include "hacc_fft.h"

class HaccFFTForwardSolve : public HaccFFT {

    public:
        
        HaccFFTForwardSolve(
                    SWMUserIF* user_if,

                    bool* done,

                    HaccConfig & config,
                    double buffer_copy_MBps,
                    double fft_work_per_second
                    );

        void call();

};

#endif
