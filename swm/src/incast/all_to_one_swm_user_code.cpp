#include "all_to_one_swm_user_code.h"

AllToOneSWMUserCode::AllToOneSWMUserCode(
    boost::property_tree::ptree cfg,
    void**& generic_ptrs
) :
    dst_rank_id(cfg.get<uint32_t>("dst_rank_id",0)),
    scattered_start(cfg.get<bool>("scattered_start", false)),
    start_delay_max(cfg.get<uint32_t>("start_delay_max", 0)),
    synchronous(cfg.get<bool>("synchronous", 0)),
    use_any_src(cfg.get<bool>("use_any_src", 0)),
    blocking_comm(cfg.get<bool>("blocking_comm", 0)),
    debug(cfg.get<bool>("debug", false))
{

    // extract the src/dst rank id intervals
    int num = 0;
    BOOST_FOREACH(const boost::property_tree::ptree::value_type &v, cfg.get_child("src_rank_id_interval"))
    {
        std::string value = v.second.data();

        if(num == 0) min_source_id = atoi(value.c_str());
        if(num == 1) max_source_id = atoi(value.c_str());

        num++;
    }
    assert(num == 2);


    assert(dst_rank_id < process_cnt);
}

void
AllToOneSWMUserCode::call()
{

  uint32_t *send_handles = NULL;
  uint32_t *recv_handles = NULL;

  uint32_t send_limit = 1;
  uint32_t recv_limit = (max_source_id - min_source_id) + 1;

  //SWMPiggybackBase* dummy_piggyback = nullptr;

  if(synchronous)
    {
      send_handles = new uint32_t[send_limit * iteration_cnt];
      recv_handles = new uint32_t[recv_limit * iteration_cnt];
    }


    if ((process_id != dst_rank_id) && (process_id >= min_source_id && process_id <= max_source_id) )   // do not send messages to self
    {

        for(uint32_t iter=0; iter < iteration_cnt; iter++)
        {

            //msg_traffic_desc msg_desc;

            //GetMsgDetails(&msg_desc);

            // if we want to scatter the start time, we mimic this delay with a compute delay
            if(scattered_start)
              {
                assert(start_delay_max > 0);
		/* TODO: Use a better random number generator here. */
                uint32_t start_delay = rand() % start_delay_max;
                std::cout << "process_id: " << process_id << " delay start by " << start_delay << " cycles" << std::endl;
                SWM_Compute(start_delay);
              }

            /*if(!synchronous)
              {

                SWM_Synthetic(
                              dst_rank_id,  //dst
                              msg_desc.msg_req_vc,
                              msg_desc.msg_rsp_vc,
                              msg_desc.pkt_rsp_vc,
                              msg_desc.msg_req_bytes,
                              msg_desc.msg_rsp_bytes,
                              msg_desc.pkt_rsp_bytes,
                              msg_desc.msg_req_routing_type,
                              msg_desc.msg_rsp_routing_type,
                              msg_desc.pkt_rsp_routing_type,
                              dummy_piggyback, //NULL,
                              msg_desc.attribute
#ifdef FABSIM_EMULATION
                              , msg_desc.l2_encoding
#endif
                              );


                if(debug)
                  {
                    std::cout << "process_id: " << process_id << " sent synthetic message to destination: " << dst_rank_id << ", iter: " << iter << " @ "  << SWM_Clock() << std::endl;
                  }

              }
            else
              {*/
                
                //uint32_t process_id_offset = ( (process_id + 1) << 32);
                //uint32_t iter_offset       = ( (iter + 1) << 8);
                //SWM_TAG this_tag = SWM_APP_TAG_BASE + process_id_offset + iter_offset;
                uint32_t iter_offset = (process_cnt * (iter) );
                SWM_TAG this_tag = SWM_APP_TAG_BASE + (sizeof(SWM_TAG) * ( (process_id + 1) + iter_offset) ); //(iter+1) );
                //uint32_t send_handle[send_limit];
                uint32_t send_count = 0;

                if(!blocking_comm)
                  {

                    SWM_Isend(
                              dst_rank_id,
                              SWM_COMM_WORLD,
                              this_tag,
			      -1, 
			      -1,
                              NO_BUFFER,
			      0, 
			      0,
                              &(send_handles[send_count]),
                              0,
                              0
                              );
                  }
                else
                  {
                    SWM_Send(
                             dst_rank_id,
                             SWM_COMM_WORLD,
                             this_tag,		
			     -1,// req-vc
			     -1, //resp-vc
                             NO_BUFFER,
                             0, //req-bytes
                             0, //resp-bytes
                             0,//routing type
                             0 //routing type
                             );
                  }

                if(!blocking_comm)
                  {
                    SWM_Waitall(send_limit, send_handles);
                  }

                if(debug)
                  {
                    std::cout << "process_id: " << process_id << " sent message to destination: " << dst_rank_id << ", tag: " << this_tag << ", iter: " << iter  << std::endl;
                  }

              //} // else(synchronous) 
	    //MM comment: no def for SWM_Noop in codes
            /*for(uint32_t noop=0; noop<noop_cnt; noop++)
            {
                SWM_Noop();
            }*/
            if (compute_delay)
                SWM_Compute(compute_delay);

        } // end-for(iteration_cnt)
    }
    else if(synchronous && (process_id == dst_rank_id) )
      {

        // need to receive from everybody every iteration...
        for(uint32_t iter = 0; iter < iteration_cnt; iter++)
          {

            uint32_t count = 0;
            
            for(uint32_t index = min_source_id; index <= max_source_id; index++, count++)
              {
                
                uint32_t iter_offset = (process_cnt * (iter) );
                //SWM_TAG this_tag = SWM_APP_TAG_BASE + (sizeof(SWM_TAG) * (index + 1) * (iter+1) );
                SWM_TAG this_tag = SWM_APP_TAG_BASE + (sizeof(SWM_TAG) * ( (index + 1) + iter_offset) );

                uint32_t receive_from_proc = (!use_any_src) ? index : -1;
                
                if(debug)
                  {
                    std::cout << "process_id: " << process_id << " expecting to recv data from: " << receive_from_proc << " with recv tag: " << this_tag << " | iter_" << iter << std::endl;
                  }
                

                if(!blocking_comm)
                  {
                    SWM_Irecv(
                              receive_from_proc,
                              SWM_COMM_WORLD,
                              this_tag,
                              NO_BUFFER,
                              &(recv_handles[count])
                              );
                  }
                else
                  {
                    SWM_Recv(
                             receive_from_proc,
                             SWM_COMM_WORLD,
                             this_tag,
                             NO_BUFFER
                             );
                  }

                if(debug)
                  {
                    std::cout << "process_id: " << process_id << " received data from src: " << index << ", iteration: " << iter  << std::endl;
                  }

              } // end of for-loop(all_sources)
        
            if(!blocking_comm)
              {
                SWM_Waitall(recv_limit, recv_handles);
              }

          } // end for-loop(iteration_cnt)

      } // end of else if(synchronous && (process_id == dst_rank_id) )

    SWM_Finalize();
    assert(0);

}

