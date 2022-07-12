/*
 * =====================================================================================
 *
 *       Filename:  layered_allbroadcast.h
 *
 *    Description:
 *
 *        Version:  1.0
 *        Created:  6/22/2021
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Neil McGlohon
 *        Company:  Rensselaer Polytechnic Institute
 *
 * =====================================================================================
 */

#ifndef _LAYERED_ALL_BROADCAST_TEMPLATE_USER_CODE_
#define _LAYERED_ALL_BROADCAST_TEMPLATE_USER_CODE_

#define SWM_APP_TAG_BASE 0

#include <boost/property_tree/ptree.hpp>

#include <string>
#include <iostream>
#include <random>
#include <algorithm>
#include <vector>
#include <regex>

#include "swm-include.h"
using namespace std;

class LayeredAllBroadcast
{

public:

    LayeredAllBroadcast(
//        SWMUserIF* user_if,
        boost::property_tree::ptree cfg,
        void**& generic_ptrs
    );

    void call();

protected:


    uint32_t process_id;
    uint32_t process_cnt;
    uint32_t iteration_cnt;

    uint32_t total_layers;
    uint32_t initial_layer_size;
    double layer_growth_rate;
    double grad_compression_rate;

    // use blocking (Send/Recv)
    bool blocking_comm;

    // for debugging
    bool show_iterations;
    bool debug;


private:

    void execute_comp_gradient_comm(int current_layer);
    void execute_weights_comm(int current_layer);
};

#endif
