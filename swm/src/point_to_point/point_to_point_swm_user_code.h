/*
 * =====================================================================================
 *
 *       Filename:  point_to_point_swm_user_code.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  07/10/2013 11:11:02 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  John Thompson (FPD), john.d.thompson@intel.com
 *        Company:  Intel
 *
 * =====================================================================================
 */

#ifndef _POINT_TO_POINT_SWM_USER_CODE_
#define _POINT_TO_POINT_SWM_USER_CODE_

#include <boost/property_tree/ptree.hpp>
#include "swm-include.h"
class PointToPointSWMUserCode {

    public:

        PointToPointSWMUserCode(
                boost::property_tree::ptree cfg,
                void**& generic_ptrs
                );

        void call(int process_id);

        bool IsMessaging(int process_id);

    protected:
        uint32_t src_rank_id;
        uint32_t dst_rank_id;
        bool two_sided;
        int iteration_cnt;
};

#endif
