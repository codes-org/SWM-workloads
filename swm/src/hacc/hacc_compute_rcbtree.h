#ifndef _HACC_COMPUTE_RCBTREE_HPP
#define _HACC_COMPUTE_RCBTREE_HPP

#include <stdio.h>

#include <boost/random/normal_distribution.hpp>
//#include <boost/random/mersenne_twister.hpp>  //BOZO -- this is causing issues

#include "hacc_config.h"

#include "swm_user_code.h"
#include "swm.h"
#include "swm_process_app_if.h"
#include "app_base_swm_user_code.h"

class HaccComputeRCBTree : public SWMUserCode {

    public:

        HaccComputeRCBTree(

                SWMUserIF* user_if,

                bool* done_from_parent,

                HaccConfig & config,
                double nint_mean,
                double nint_delta,
                double nint_per_wall_second
                );

    void call();

    protected:

        bool* done_to_parent;

        HaccConfig & config;

        // Ensemble mean (across ranks) of number of force interactions/rank
        double nint_mean;

        // Standard deviation (across ranks) of number of force interactions/rank,
        // normalized by nint_mean (i.e. delta=sigma/mu)
        double nint_delta;

        // Number of interactions computed per wall second
        double nint_per_wall_second;

        double nint;


};

#endif
