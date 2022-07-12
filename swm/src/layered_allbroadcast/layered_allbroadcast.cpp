#include "layered_allbroadcast.h"
#include "math.h"
#include <stdarg.h> //log printf wrapper

bool stdout_log = false;

static void print_log(const char *format, ...)
{
    va_list args;
    va_start(args, format);

    if(stdout_log)
        vprintf(format, args);

    va_end(args);
}


LayeredAllBroadcast::LayeredAllBroadcast(
		boost::property_tree::ptree cfg,
		void**& generic_ptrs
		) :
	process_cnt(cfg.get<uint32_t>("jobs.size", 1)),
	iteration_cnt(cfg.get<uint32_t>("jobs.cfg.iteration_cnt", 1)),
	total_layers(cfg.get<uint32_t>("jobs.cfg.total_layers",50)),
	initial_layer_size(cfg.get<double>("jobs.cfg.initial_layer_size",8192)),
	layer_growth_rate(cfg.get<double>("jobs.cfg.layer_growth_rate",1.12)),
	grad_compression_rate(cfg.get<double>("jobs.cfg.first_comm_compression_rate",32)),
	blocking_comm(cfg.get<bool>("jobs.cfg.blocking_comm", false)),
	show_iterations(cfg.get<bool>("jobs.cfg.show_iterations", false)),
	debug(cfg.get<bool>("jobs.cfg.debug", false))
{
	stdout_log = debug;
	process_id = *((int*)generic_ptrs[0]);
}

void LayeredAllBroadcast::call()
{

	int iter_marker = 0;
	for (int iter = 0; iter < iteration_cnt; iter++)
	{
		if (show_iterations) {
			SWM_Mark_Iteration(iter_marker);
			iter_marker++;
		}

		for (int i = 0; i < total_layers; i++)
		{
			if (process_id == 0) print_log("LayeredAllBcast Layer %d Comp Grad\n",i);
			execute_comp_gradient_comm(i);
			if (process_id == 0) print_log("LayeredAllBcast Layer %d Weights\n",i);
			execute_weights_comm(i);
		}
		// SWM_Allreduce(32, 0, SWM_COMM_WORLD, -1, -1, NO_BUFFER, NO_BUFFER);
		// SWM_Barrier(SWM_COMM_WORLD, -1, -1, NO_BUFFER, 0, 0, 0, 0);
		if (process_id == 0) print_log("LayeredAllBcast Iteration %d/%d Completed\n",iter+1,iteration_cnt);

		if (show_iterations) {
			SWM_Mark_Iteration(iter_marker);
			iter_marker++;
		}
	}


	SWM_Finalize();
}


void LayeredAllBroadcast::execute_comp_gradient_comm(int current_layer)
{
	double grad_size = (initial_layer_size * (pow(layer_growth_rate,current_layer)))/grad_compression_rate;
	double piece_size = grad_size / process_cnt;

    uint32_t *h = (uint32_t*)calloc(process_cnt-1, sizeof(uint32_t));
	uint32_t *h2 = (uint32_t*)calloc(process_cnt-1, sizeof(uint32_t));
	int send_count = 0;
	int recv_count = 0;

	for(int i = 0; i < process_cnt; i++)
	{
		if (i != process_id)
		{
			SWM_Irecv(i,SWM_COMM_WORLD, 0, NO_BUFFER, &(h2[recv_count]));
			recv_count++;
		}
	}

	for(int i = 0; i < process_cnt; i++)
	{
		if (i != process_id)
		{
			SWM_Isend(i,SWM_COMM_WORLD, 0, -1, -1, NO_BUFFER, (int)piece_size, 0, &(h[send_count]),0,0);
			send_count++;
		}
	}



	SWM_Waitall(send_count, h);
	SWM_Waitall(recv_count, h2);
	free(h);
	free(h2);
}

void LayeredAllBroadcast::execute_weights_comm(int current_layer)
{
	double weights_size = (initial_layer_size * (pow(layer_growth_rate,current_layer)));
	double piece_size = weights_size / process_cnt;

    uint32_t *h = (uint32_t*)calloc(process_cnt-1, sizeof(uint32_t));
	uint32_t *h2 = (uint32_t*)calloc(process_cnt-1, sizeof(uint32_t));
	int send_count = 0;
	int recv_count = 0;
	for(int i = 0; i < process_cnt; i++)
	{
		if (i != process_id)
		{
			SWM_Irecv(i,SWM_COMM_WORLD, 1, NO_BUFFER, &(h2[recv_count]));
			recv_count++;
		}
	}

	for(int i = 0; i < process_cnt; i++)
	{
		if (i != process_id)
		{
			SWM_Isend(i,SWM_COMM_WORLD, 1, -1, -1, NO_BUFFER, (int)piece_size, 0, &(h[send_count]),0,0);
			send_count++;
		}
	}

	SWM_Waitall(send_count, h);
	SWM_Waitall(recv_count, h2);
	free(h);
	free(h2);
}



/*
 * Local variables:
 *  c-indent-level: 4
 *  c-basic-offset: 4
 * End:
 *
 * vim: ft=c ts=8 sts=4 sw=4 expandtab
 */
