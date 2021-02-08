/*
 * =====================================================================================
 *
 *       Filename:  spread_swm_user_code.h
 *
 *    Description:
 *
 *        Version:  1.0
 *        Created:  09/26/2020 01:05:02 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Kevin A. Brown, kb@anl.gov
 *        Company:  Argonne Nat Lab
 *
 * =====================================================================================
 */

#ifndef _ALLREDUCE_TEMPLATE_USER_CODE_
#define _ALLREDUCE_TEMPLATE_USER_CODE_

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

class AllReduceSWMUserCode 
{

public:

    AllReduceSWMUserCode(
//        SWMUserIF* user_if,
        boost::property_tree::ptree cfg,
        void**& generic_ptrs
    );

    void call();

protected:
    uint32_t request_vc;
    uint32_t response_vc;
    uint32_t msg_req_bytes;
    uint32_t msg_rsp_bytes;

    uint32_t process_id;
    uint32_t process_cnt;
    uint32_t iteration_cnt;
    uint32_t compute_delay;

    // for debugging
    bool show_iterations;
    bool debug;

};

#endif
