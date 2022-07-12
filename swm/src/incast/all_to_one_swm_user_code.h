/*
 * =====================================================================================
 *
 *       Filename:  all_to_one_swm_user_code.h
 *
 *    Description:
 *
 *        Version:  1.0
 *        Created:  12/3/2013 01:05:02 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Nate Andrysco, nathan.r.andrysco@intel.com
 *        Company:  Intel
 *
 * =====================================================================================
 */

#ifndef _ALL_TO_ONE_TEMPLATE_USER_CODE_
#define _ALL_TO_ONE_TEMPLATE_USER_CODE_

#define SWM_APP_TAG_BASE 0

#include <boost/property_tree/ptree.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/foreach.hpp>

#include <string>
#include <iostream>
#include <random>
#include <algorithm>
#include <vector>
#include <regex>

#include "swm-include.h"
using namespace std;

class AllToOneSWMUserCode 
{

public:

    AllToOneSWMUserCode(
//        SWMUserIF* user_if,
        boost::property_tree::ptree cfg,
        void**& generic_ptrs
    );

    void call();

protected:
    std::string req_vcs_string;
    std::string rsp_vcs_string;
    uint32_t msg_req_bytes;
    uint32_t msg_rsp_bytes;
    uint32_t pkt_rsp_bytes;
    uint32_t dst_rank_id;
    uint32_t msg_size; // MM addition


    uint32_t process_id;
    uint32_t process_cnt;
    uint32_t iteration_cnt;
    uint32_t noop_cnt;
    uint32_t compute_delay;

    std::vector<uint32_t> req_vcs;
    std::vector<uint32_t> rsp_vcs;

    uint32_t min_source_id;
    uint32_t max_source_id;

    bool randomize_comm_order;

    // are we staggering the start time of the srcs
    bool scattered_start;

    // if using staggered start delay, this is the max used in the RNG when computing delay
    uint32_t start_delay_max;

    // use isend/irecv instead of synthetic                                                                                                                                   
    bool synchronous;

    // use __ANY__ at the receive side for synchronous                                                                                                                
    bool use_any_src;

    // use blocking (Send/Recv)                                                                                                                              
    bool blocking_comm;

    // for debugging
    bool show_iterations;
    bool debug;

};

#endif
