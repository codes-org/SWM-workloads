#include "allreduce.h"

AllReduceSWMUserCode::AllReduceSWMUserCode(
		boost::property_tree::ptree cfg,
		void**& generic_ptrs
		) :
	process_cnt(cfg.get<uint32_t>("jobs.size", 1)),
	iteration_cnt(cfg.get<uint32_t>("jobs.cfg.iteration_cnt", 1)),
	msg_req_bytes(cfg.get<uint32_t>("jobs.cfg.msg_req_bytes", 1024)),
	msg_rsp_bytes(cfg.get<uint32_t>("jobs.cfg.msg_rsp_bytes", 0)),
	compute_delay(cfg.get<uint32_t>("jobs.cfg.compute_delay", 0)),
        show_iterations(cfg.get<bool>("jobs.cfg.show_iterations", false))
{

	request_vc = 0;
	response_vc = 0;

	process_id = *((int*)generic_ptrs[0]);
}

void
AllReduceSWMUserCode::call()
{
	/* Print job description */
	if(process_id == 0)
	{
		std::cout << std::endl << "JOB: Allreduce | size: " << process_cnt;
		std::cout << " | interation_cnt: " << iteration_cnt;
		std::cout << " | compute_delay: " << compute_delay << std::endl;
	}


	uint32_t tag = 0;
	for(uint32_t iter=0; iter < iteration_cnt; iter++)
	{

		if (compute_delay)
			SWM_Compute(compute_delay);

		//if(process_id == 0)
		//{
		    /* Print the start time of the Allreduce on the rank */
                    if(show_iterations){
		        SWM_Mark_Iteration(tag);
			tag = tag +1;
                    }
		//}

		SWM_Allreduce(
				msg_req_bytes, // payload
				msg_rsp_bytes, // pkt_rsp_bytes
				SWM_COMM_WORLD, 
				request_vc,
				response_vc,
				NO_BUFFER,
				NO_BUFFER);

		//if(process_id == 0)
		//{
		    /* Print the end time of the Allreduce call on the rank */
                    if(show_iterations){
			SWM_Mark_Iteration(tag);        
			tag = tag +1;
                    }
		//}
	}



	SWM_Finalize();
}

/*
 * Local variables:
 *  c-indent-level: 4
 *  c-basic-offset: 4
 * End:
 *
 * vim: ft=c ts=8 sts=4 sw=4 expandtab
 */
