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

#include "app_base_swm_user_code.h"

class PointToPointSWMUserCode : public AppBaseSWMUserCode {

    public:

        PointToPointSWMUserCode(
                SWMUserIF* user_if,
                boost::property_tree::ptree cfg,
                void**& generic_ptrs
                );

        void call();

        bool IsMessaging();

    protected:
        uint32_t src_rank_id;
        uint32_t dst_rank_id;
        bool two_sided;

};

#endif
