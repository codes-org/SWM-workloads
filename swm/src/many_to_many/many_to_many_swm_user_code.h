/*
 * =====================================================================================
 *
 *       Filename:  many_to_many_swm_user_code.h
 *
 *    Description:
 *
 *         Author:  Kevin Brown, kabrown@anl.gov
 *
 * =====================================================================================
 */

#ifndef _MANY_TO_MANY_TEMPLATE_USER_CODE_
#define _MANY_TO_MANY_TEMPLATE_USER_CODE_

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

class ManyToManySWMUserCode 
{

public:

    ManyToManySWMUserCode(
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

    uint32_t process_id;
    uint32_t process_cnt;
    uint32_t iteration_cnt;
    uint32_t noop_cnt;
    uint32_t compute_delay;

    std::vector<uint32_t> req_vcs;
    std::vector<uint32_t> rsp_vcs;

    uint32_t min_src_id;
    uint32_t max_src_id;
    uint32_t min_dst_id;
    uint32_t max_dst_id;

    bool randomize_comm_order;

    // are we staggering the start time of the srcs
    bool scattered_start;

    // one-to-one pairs between sender and recv?
    bool fixed_pairs;

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
