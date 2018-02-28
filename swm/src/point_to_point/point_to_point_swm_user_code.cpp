#include "point_to_point_swm_user_code.h"

PointToPointSWMUserCode::PointToPointSWMUserCode(
        boost::property_tree::ptree cfg,
        void**& generic_ptrs
        ) :
    src_rank_id(cfg.get<uint32_t>("src_rank_id",0)),
    dst_rank_id(cfg.get<uint32_t>("dst_rank_id",0)),
    two_sided(cfg.get<bool>("two_sided",0))
{ 
}

bool
PointToPointSWMUserCode::IsMessaging(int process_id) {
    return (process_id == src_rank_id);
}

void
PointToPointSWMUserCode::call(int process_id) 
{

    SWM_Init();

    if(process_id == src_rank_id)
    {
        for(uint32_t iter=0; iter<iteration_cnt; iter++) 
        {

            msg_traffic_desc msg_desc;

            GetMsgDetails(&msg_desc);

            if((iter + 1) == iteration_cnt) {
                std::cout << "set sent_last_message" << std::endl;
                sent_last_message=true;
            }

            SWM_Synthetic(
                    dst_rank_id,
                    msg_desc.msg_req_vc,
                    msg_desc.msg_rsp_vc,
                    msg_desc.pkt_rsp_vc,
                    msg_desc.msg_req_bytes,
                    msg_desc.msg_rsp_bytes,
                    msg_desc.pkt_rsp_bytes,
                    msg_desc.msg_req_routing_type,
                    msg_desc.msg_rsp_routing_type,
                    msg_desc.pkt_rsp_routing_type,
                    NULL,
                    msg_desc.attribute
#ifdef FABSIM_EMULATION
                    , msg_desc.l2_encoding
#endif
                    );
            for(uint32_t noop=0; noop<noop_cnt; noop++) {
                SWM_Noop();
            }
            if (compute_delay)
                SWM_Compute(compute_delay);
        }
    }
    else if(two_sided && (process_id == dst_rank_id))
    {
        for(uint32_t iter=0; iter<iteration_cnt; iter++) 
        {

            msg_traffic_desc msg_desc;

            GetMsgDetails(&msg_desc);

            SWM_Synthetic(
                    src_rank_id,
                    msg_desc.msg_req_vc,
                    msg_desc.msg_rsp_vc,
                    msg_desc.pkt_rsp_vc,
                    msg_desc.msg_req_bytes,
                    msg_desc.msg_rsp_bytes,
                    msg_desc.pkt_rsp_bytes,
                    msg_desc.msg_req_routing_type,
                    msg_desc.msg_rsp_routing_type,
                    msg_desc.pkt_rsp_routing_type,
                    NULL,
                    msg_desc.attribute
#ifdef FABSIM_EMULATION
                    , msg_desc.l2_encoding
#endif
                    );
            for(uint32_t noop=0; noop<noop_cnt; noop++) {
                SWM_Noop();
            }
            if (compute_delay)
                SWM_Compute(compute_delay);
        }
    }

    SWM_Finalize();
    assert(0);

}

