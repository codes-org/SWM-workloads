#include "many_to_many_swm_user_code.h"

ManyToManySWMUserCode::ManyToManySWMUserCode(
		boost::property_tree::ptree cfg,
		void**& generic_ptrs
		) :
	process_cnt(cfg.get<uint32_t>("jobs.size", 1)),
	iteration_cnt(cfg.get<uint32_t>("jobs.cfg.iteration_cnt", 1)),
	msg_req_bytes(cfg.get<uint32_t>("jobs.cfg.msg_req_bytes", 0)),
	msg_rsp_bytes(cfg.get<uint32_t>("jobs.cfg.msg_rsp_bytes", 0)),
	compute_delay(cfg.get<uint32_t>("jobs.cfg.compute_delay", 0)),
	use_any_src(cfg.get<bool>("jobs.cfg.use_any_src", false)),
	blocking_comm(cfg.get<bool>("jobs.cfg.blocking_comm", false)),
	scattered_start(cfg.get<bool>("jobs.cfg.scattered_start", false)),
	fixed_pairs(cfg.get<bool>("jobs.cfg.fixed_pairs", false)),
	start_delay_max(cfg.get<uint32_t>("jobs.cfg.start_delay_max", 0)),
	randomize_comm_order(cfg.get<bool>("jobs.cfg.randomize_communication_order", false)),
        show_iterations(cfg.get<bool>("jobs.cfg.show_iterations", false)),
	debug(cfg.get<bool>("jobs.cfg.debug", false))
{

	// extract the src/dst rank id intervals
	int num = 0;
	BOOST_FOREACH(const boost::property_tree::ptree::value_type &v, cfg.get_child("jobs.cfg.src_rank_id_interval"))
	{
		std::string value = v.second.data();

		if(num == 0) min_src_id = atoi(value.c_str());
		if(num == 1) max_src_id = atoi(value.c_str());

		num++;
	}
	assert(num == 2);

	// extract the src/dst rank id intervals
	num = 0;
	BOOST_FOREACH(const boost::property_tree::ptree::value_type &v, cfg.get_child("jobs.cfg.dst_rank_id_interval"))
	{
		std::string value = v.second.data();

		if(num == 0) min_dst_id = atoi(value.c_str());
		if(num == 1) max_dst_id = atoi(value.c_str());

		num++;
	}
	assert(num == 2);

	assert(max_src_id < process_cnt);
	assert(max_dst_id < process_cnt);

	process_id = *((int*)generic_ptrs[0]);
}

void
ManyToManySWMUserCode::call()
{
	/* Print job description */
	if(process_id == 0)
	{
		std::cout << std::endl << "JOB: Bulk_data | size: " << process_cnt;
		std::cout << " | interation_cnt: " << iteration_cnt;
		std::cout << " | msg_req_bytes: " << msg_req_bytes;
		std::cout << " | msg_rsp_bytes: " << msg_rsp_bytes;
		std::cout << " | src_rank_id_interval: " << min_src_id << "-" << max_src_id;
		std::cout << " | dst_rank_id_interval: " << min_dst_id << "-" << max_dst_id;
		std::cout << " | scattered_start: " << scattered_start;
		std::cout << " | compute_delay: " << compute_delay << std::endl;
	}
	uint32_t *send_handles = NULL;
	uint32_t *recv_handles = NULL;

	uint32_t send_limit = (max_dst_id - min_dst_id) + 1;
	uint32_t recv_limit = (max_src_id - min_src_id) + 1;

	if (fixed_pairs)
	{
		send_limit = 1;
		recv_limit = 1;
	}

	send_handles = new uint32_t[send_limit * iteration_cnt];
	recv_handles = new uint32_t[recv_limit * iteration_cnt];


	if (process_id >= min_src_id && process_id <= max_src_id) // Sending processes
	{
		// if we want to scatter the start time, we mimic this delay with a compute delay
		if(scattered_start)
		{
			assert(start_delay_max > 0);
			/* TODO: Use a better random number generator here. */
			uint32_t start_delay = rand() % start_delay_max;
			std::cout << std::endl << "process_id: " << process_id << " delay start by " << start_delay << " cycles";
			SWM_Compute(start_delay);
		}

		uint32_t marker = 0;
		for(uint32_t iter=0; iter < iteration_cnt; iter++)
		{
			if (compute_delay)
				SWM_Compute(compute_delay);


			//if(show_iterations)
                        //  SWM_Mark_Iteration(marker);
			//  marker++;
                        //}
			uint32_t send_count = 0;

			uint32_t curr_target;
			if (fixed_pairs)
			{
				curr_target = (process_id - min_src_id) + min_dst_id;
			}
			else
			{
				curr_target = (process_id % send_limit) + min_dst_id;
			}
			if (curr_target > max_dst_id)
			{
				std::cout << std::endl << "process_id: " << process_id << " - unused source. stopping.";
				return;
			}

			for(uint32_t sent=0; sent < send_limit; sent++, send_count++)
			{
				//uint32_t process_id_offset = ( (process_id + 1) << 32);
				//uint32_t iter_offset       = ( (iter + 1) << 8);
				//SWM_TAG this_tag = SWM_APP_TAG_BASE + process_id_offset + iter_offset;
				//uint32_t iter_offset = (process_cnt * (iter) );
				//SWM_TAG this_tag = SWM_APP_TAG_BASE + (sizeof(SWM_TAG) * ( (process_id + 1) + iter_offset) ); //(iter+1) );
				SWM_TAG this_tag = SWM_APP_TAG_BASE;
				//uint32_t send_handle[send_limit];

				SWM_Isend(
						curr_target,
						SWM_COMM_WORLD,
						this_tag,
						-1, 
						-1,
						NO_BUFFER,
						msg_req_bytes, 
						msg_rsp_bytes,
						&(send_handles[send_count]),
						0,
						0
					 );
				if(debug)
				{
					std::cout << std::endl << "process_id: " << process_id << " sent message to destination: " << curr_target << ", tag: " << this_tag << ", iter: " << iter ;
				}

				if (!fixed_pairs)
				{
					curr_target++;
					if (curr_target > max_dst_id)
					{
						curr_target = min_dst_id;
					}
				}
			}

			SWM_Waitall(send_limit, send_handles);


                        if(show_iterations){
			    SWM_Mark_Iteration(marker);
			    marker++;
                        }
		} // end-for(iteration_cnt)
	}
	else if (process_id >= min_dst_id && process_id <= max_dst_id) // Recieving processes
	{
		uint32_t     receive_from_proc =  -1;

		// Are were using fixed sender-reciever pairs? Or do each reciever get data from all senders?
		if (fixed_pairs)
		{
			min_src_id = (process_id - min_dst_id) + min_src_id;
			if(min_src_id > max_src_id){
				std::cout << std::endl << "process_id: " << process_id << " - unused destination. stopping.";
				return;
			}
			max_src_id = min_src_id;
		}
		// need to receive from each sender every iteration...
		for(uint32_t iter = 0; iter < iteration_cnt; iter++)
		{
			uint32_t count = 0;

			for(uint32_t index = min_src_id; index <= max_src_id; index++, count++)
			{
				//uint32_t iter_offset = (process_cnt * (iter) );
				//SWM_TAG this_tag = SWM_APP_TAG_BASE + (sizeof(SWM_TAG) * (index + 1) * (iter+1) );
				//SWM_TAG this_tag = SWM_APP_TAG_BASE + (sizeof(SWM_TAG) * ( (index + 1) + iter_offset) );
				SWM_TAG this_tag = SWM_APP_TAG_BASE;

				//uint32_t receive_from_proc = (!use_any_src) ? index : -1;

				if(debug)
				{
					std::cout  << std::endl << "process_id: " << process_id << " expecting to recv data from: " << index  << " with recv tag: " << this_tag << " | iter_" << iter;
				}


				SWM_Irecv(
						index,
						SWM_COMM_WORLD,
						this_tag,
						NO_BUFFER,
						&(recv_handles[count])
					 );

			} // end of for-loop(all_sources)

			SWM_Waitall(recv_limit, recv_handles);
			if(debug)
			{
				std::cout << std::endl << "process_id: " << process_id << " received all data. iteration: " << iter ;
			}

			//if(show_iterations)
                        //  SWM_Mark_Iteration(iter);
		} // end for-loop(iteration_cnt)

	} // end of else if(synchronous && (process_id == dst_rank_id) )

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
