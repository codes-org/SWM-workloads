
#include "hacc_compute_rcbtree.h"

HaccComputeRCBTree::HaccComputeRCBTree(

        SWMUserIF* user_if,

        bool* done_from_parent,

        HaccConfig & config,
        double nint_mean,
        double nint_delta,
        double nint_per_wall_second
        ) :
    SWMUserCode(user_if),
    done_to_parent(done_from_parent),
    config(config),
    nint_mean(nint_mean),
    nint_delta(nint_delta),
    nint_per_wall_second(nint_per_wall_second)
{

    if (nint_delta > 0.0) {
        // Draw a Gaussian random sample using boost::random
        //double nint_sigma = nint_mean*nint_delta;

        // Make sure to seed the RNG with our rank id, so that a) it's
        // reproducible and b) all the ranks don't get the same variate
        // BOZO -- this is causing trouble
        //boost::mt19937 rng(config.myrank);
        //boost::normal_distribution<double> d(nint_mean, nint_sigma);
        //nint = d(rng);
        nint = nint_mean;
    } else {
        // Just use the mean
        nint = nint_mean;
    }
}

void
HaccComputeRCBTree::call() { //build_tree_and_evaluate_forces() {
    if (enable_contexts)
    while(1) {
        // TODO:compute: tree build

        // Compute interactions
        //backend.compute_seconds(nint/nint_per_wall_second);
        *done_to_parent = false;
        SWM_Compute(nint/nint_per_wall_second);
        *done_to_parent = true; yield();
    }
    else
    {
        *done_to_parent = false;
        SWM_Compute(nint/nint_per_wall_second);
        *done_to_parent = true; yield();
    }
}
