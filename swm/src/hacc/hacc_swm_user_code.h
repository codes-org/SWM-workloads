#ifndef _HACC_TEMPLATE_USER_CODE_
#define _HACC_TEMPLATE_USER_CODE_

#include <boost/property_tree/json_parser.hpp>
#include "app_base_swm_user_code.h"

#include "defines.h"
#include "resolve_string_environment_variables.h"


#include "hacc_config.h"
#include "hacc_compute_rcbtree.h"
#include "hacc_fft_forward_solve.h"
#include "hacc_fft_backward_solve_gradient.h"
#include "hacc_exchange.h"
#include "hacc_timestep.h"

class HaccTimestep;

class HACCSWMUserCode : public SWMUserCode {


    public:

        HACCSWMUserCode(
                SWMUserIF* user_if,
                boost::property_tree::ptree cfg,
                void**& generic_ptrs
                );

        ~HACCSWMUserCode() {
            if (timestep)
                delete timestep;
        }

        void call();

    protected:

        SWMUserIF* user_if;

        uint32_t request_vc;
        uint32_t response_vc;
        uint32_t pkt_rsp_bytes;

        std::string gen_cfg_filename;
        boost::property_tree::ptree cfg;
        boost::property_tree::ptree gen_cfg;

        bool done_to_child;

        int ng;
        int nranks; //8;

        int rank_shape_3d  [3];
        int rank_shape_2d_x[3];
        int rank_shape_2d_y[3];
        int rank_shape_2d_z[3];

        double box_length;

        HaccTimestep* timestep = nullptr;
};

#endif
