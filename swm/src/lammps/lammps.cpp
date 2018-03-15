#include <time.h>
#include <math.h>
#include <inttypes.h>

#include "lammps.h"

static double round(double x, int p)
{
    return floor(x*pow(10,p) + 0.5)/pow(10,p);
}

/*MM: Function signature modified SWMUserIF removed. */
LAMMPS_SWM::LAMMPS_SWM(
    boost::property_tree::ptree cfg,
    void**& generic_ptrs
) :
    process_cnt(cfg.get<uint32_t>("jobs.size", 1)),
    x_rep(cfg.get<uint32_t>("jobs.cfg.num_x_replicas", 1)),
    y_rep(cfg.get<uint32_t>("jobs.cfg.num_y_replicas", 1)),
    z_rep(cfg.get<uint32_t>("jobs.cfg.num_z_replicas", 1)),
    num_timesteps(cfg.get<uint32_t>("jobs.cfg.num_time_steps", 100)),
    router_freq(cfg.get<double>("jobs.cfg.router_freq", 800e6)),
    cpu_freq(cfg.get<double>("jobs.cfg.cpu_freq", 1.2e9)),
    cpu_sim_speedup(cfg.get<double>("jobs.cfg.cpu_sim_speedup", 1.0))
{

    /*TODO MM: Figure out if we need to keep this code. 
     * msg_traffic_desc msg_desc;
    GetMsgDetails(&msg_desc);

    req_vc  = msg_desc.msg_req_vc;
    resp_vc = msg_desc.msg_rsp_vc;
    rsp_bytes   = msg_desc.pkt_rsp_bytes;
    */

    /* TODO MM changes. VCs are not required for the workload part. */
    req_vc = 0;
    resp_vc = 1;
    rsp_bytes = 0;
    process_id = *((int*)generic_ptrs[0]);

    int i = 0;

    for( i = 0; i < NUM_TRANSPOSE; i++)
    {
        k_r_targets[i] = nullptr;
        k_s_targets[i] = nullptr;
        k_s_sizes[i] = nullptr;
        k_cyc[i] = 0;
        k_len[i] = 0;
    }

    gh_fw_len = gh_rw_len = k_pre_len = k_post_len = fix_len = neigh_e_len = neigh_b_len = 0;
    gh_fw_r_targets   = gh_fw_s_targets   = gh_fw_s_sizes   = nullptr;
    gh_rw_r_targets   = gh_rw_s_targets   = gh_rw_s_sizes   = nullptr;
    k_pre_r_targets   = k_pre_s_targets   = k_pre_s_sizes   = nullptr;
    k_post_r_targets  = k_post_s_targets  = k_post_s_sizes  = nullptr;
    fix_r_targets     = fix_s_targets     = fix_s_sizes     = nullptr;
    neigh_e_r_targets = neigh_e_s_targets = neigh_e_s_sizes = nullptr;
    neigh_b_r_targets = neigh_b_s_targets = neigh_b_s_sizes = nullptr;
    gh_fw_cyc   = nullptr;
    gh_rw_cyc   = nullptr;
    k_pre_cyc   = nullptr;
    k_post_cyc  = nullptr;
    fix_cyc     = nullptr;
    neigh_e_cyc = nullptr;
    neigh_b_cyc = nullptr;


    // size of the problem in (x,y,z)
    prd[0] = x_rep * (XHI_BASE - XLO_BASE);
    prd[1] = y_rep * (YHI_BASE - YLO_BASE);
    prd[2] = z_rep * (ZHI_BASE - ZLO_BASE);

    // Decompose domain to processes, store result in procNums[3]
    proc_decomposition(process_cnt, prd, procNums);

    pppm_decomposition(N_ATOMS_BASE*x_rep*y_rep*z_rep, prd, pppmGrid);

}


LAMMPS_SWM::~LAMMPS_SWM()
{
    int i = 0;

    for( i = 0; i < NUM_TRANSPOSE; i++)
    {
        if(k_r_targets[i]) delete k_r_targets[i];
        if(k_s_targets[i]) delete k_s_targets[i];
        if(k_s_sizes[i]) delete k_s_sizes[i];
    }

    if(gh_fw_r_targets) delete gh_fw_r_targets;
    if(gh_fw_s_targets) delete gh_fw_s_targets;
    if(gh_fw_s_sizes) delete gh_fw_s_sizes;
    if(gh_fw_cyc) delete gh_fw_cyc;
    if(gh_rw_r_targets) delete gh_rw_r_targets;
    if(gh_rw_s_targets) delete gh_rw_s_targets;
    if(gh_rw_s_sizes) delete gh_rw_s_sizes;
    if(gh_rw_cyc) delete gh_rw_cyc;
    if(k_pre_r_targets) delete k_pre_r_targets;
    if(k_pre_s_targets) delete k_pre_s_targets;
    if(k_pre_s_sizes) delete k_pre_s_sizes;
    if(k_pre_cyc) delete k_pre_cyc;
    if(k_post_r_targets) delete k_post_r_targets;
    if(k_post_s_targets) delete k_post_s_targets;
    if(k_post_s_sizes) delete k_post_s_sizes;
    if(k_post_cyc) delete k_post_cyc;
    if(fix_r_targets) delete fix_r_targets;
    if(fix_s_targets) delete fix_s_targets;
    if(fix_s_sizes) delete fix_s_sizes;
    if(fix_cyc) delete fix_cyc;
    if(neigh_e_r_targets) delete neigh_e_r_targets;
    if(neigh_e_s_targets) delete neigh_e_s_targets;
    if(neigh_e_s_sizes) delete neigh_e_s_sizes;
    if(neigh_e_cyc) delete neigh_e_cyc;
    if(neigh_b_r_targets) delete neigh_b_r_targets;
    if(neigh_b_s_targets) delete neigh_b_s_targets;
    if(neigh_b_s_sizes) delete neigh_b_s_sizes;
    if(neigh_b_cyc) delete neigh_b_cyc;
}

void
LAMMPS_SWM::doP2P(int len, int *r_targets, int *s_targets, int *s_sizes, long *cyc_cnt)
{
    int i = 0;
    uint32_t h = 0;
    for(i = 0; i < len; i++)
    {
        SWM_Compute(cyc_cnt[i]);
        SWM_Irecv(r_targets[i], SWM_COMM_WORLD, 0, NO_BUFFER, &h);
        SWM_Send(s_targets[i], SWM_COMM_WORLD, 0, req_vc, resp_vc, NO_BUFFER, s_sizes[i]);
        //printf("S @ %d: i %d: %d -> %d\n", process_id, i, process_id, s_targets[i]);
        //printf("W @ %d: i %d: %d -> %d\n", process_id, i, r_targets[i], process_id);
        SWM_Wait(h);
        //printf("D @ %d: i %d: %d -> %d\n", process_id, i, r_targets[i], process_id);
    }
}

void
LAMMPS_SWM::doNeighExch()
{
    int i = 0;
    uint32_t h = 0;

    // neighbor exchange
    while(i < neigh_e_len)
    {
        SWM_Compute(neigh_e_cyc[i]);
        SWM_Sendrecv(SWM_COMM_WORLD, neigh_e_r_targets[i], 0, req_vc, resp_vc, NO_BUFFER, 4, rsp_bytes, neigh_e_s_targets[i], 0, NO_BUFFER);
        if(neigh_e_r_targets[i] != neigh_e_s_targets[i])
        {
            SWM_Sendrecv(SWM_COMM_WORLD, neigh_e_s_targets[i], 0, req_vc, resp_vc, NO_BUFFER, 4, rsp_bytes, neigh_e_r_targets[i], 0, NO_BUFFER);
        }
        SWM_Irecv(neigh_e_r_targets[i], SWM_COMM_WORLD, 0, NO_BUFFER, &h);
        SWM_Send(neigh_e_s_targets[i], SWM_COMM_WORLD, 0, req_vc, resp_vc, NO_BUFFER, neigh_e_s_sizes[i]);
        SWM_Wait(h);
        i++;
        if((i < neigh_e_len) && (neigh_e_r_targets[i-1] != neigh_e_s_targets[i-1]))
        {
            SWM_Irecv(neigh_e_r_targets[i], SWM_COMM_WORLD, 0, NO_BUFFER, &h);
            SWM_Send(neigh_e_s_targets[i], SWM_COMM_WORLD, 0, req_vc, resp_vc, NO_BUFFER, neigh_e_s_sizes[i]);
            SWM_Wait(h);
            i++;
        }
    }

    // neighbor borders
    for(i = 0; i < neigh_b_len; i++)
    {
        SWM_Compute(neigh_b_cyc[i]);
        SWM_Sendrecv(SWM_COMM_WORLD, neigh_b_r_targets[i], 0, req_vc, resp_vc, NO_BUFFER, 4, rsp_bytes, neigh_b_s_targets[i], 0, NO_BUFFER);
        SWM_Irecv(neigh_b_r_targets[i], SWM_COMM_WORLD, 0, NO_BUFFER, &h);
        SWM_Send(neigh_b_s_targets[i], SWM_COMM_WORLD, 0, req_vc, resp_vc, NO_BUFFER, neigh_b_s_sizes[i]);
        SWM_Wait(h);
    }

    // allreduces
    for(i = 0; i < NUM_NEIGH_ALLREDUCE; i++)
    {
        SWM_Compute(neigh_end_cyc[i]);
        SWM_Allreduce(4, rsp_bytes, SWM_COMM_WORLD, req_vc, resp_vc, NO_BUFFER, NO_BUFFER);
    }
}

void
LAMMPS_SWM::doFFT()
{
    uint32_t *h;
    int i = 0, idx = 0;

    for(idx = 0; idx < NUM_TRANSPOSE; idx++)
    {

        h = new uint32_t[k_len[idx]];

        SWM_Compute(k_cyc[idx]);
        for(i = 0; i < k_len[idx]; i++)
        {
            SWM_Irecv(k_r_targets[idx][i], SWM_COMM_WORLD, 0, NO_BUFFER, &h[i]);
        }
        for(i = 0; i < k_len[idx]; i++)
        {
            SWM_Send(k_s_targets[idx][i], SWM_COMM_WORLD, 0, req_vc, resp_vc, NO_BUFFER, k_s_sizes[idx][i]);
        }
        SWM_Waitall(k_len[idx], h);

        delete h;
    }
}

bool
LAMMPS_SWM::neigh_check()
{
    if(neigh_check_count < NEIGH_DELAY)
    {
        neigh_check_count++;
        return false;
    }
    else
    {

        if( (neigh_check_count - NEIGH_DELAY) % NEIGH_EVERY )
        {
            neigh_check_count++;
            return false;
        }
        else
        {
            SWM_Compute(neigh_check_cyc);
            SWM_Allreduce(4, rsp_bytes, SWM_COMM_WORLD, req_vc, resp_vc, NO_BUFFER, NO_BUFFER);

            neigh_check_cumulative += neigh_check_average;

            if(neigh_check_cumulative > 1.0)
            {
                neigh_check_cumulative -= 1.0;
                neigh_check_count = 0;
                return true;
            }
        }
    }

    neigh_check_count++;
    return false;
}

/* TODO (MM): Eliminate process id from the call method. */
void
LAMMPS_SWM::call()
{
    unsigned int ts = 0;
    modelInit();

    for(ts = 0; ts < num_timesteps; ts++)
    {
        // initial integration
        SWM_Compute(start_cyc);
        SWM_Allreduce(48, rsp_bytes, SWM_COMM_WORLD, req_vc, resp_vc, NO_BUFFER, NO_BUFFER); // temperature
        SWM_Allreduce(48, rsp_bytes, SWM_COMM_WORLD, req_vc, resp_vc, NO_BUFFER, NO_BUFFER); // pressure

        // check if neighbors need to be exchanged
        if(neigh_check())
        {
            // do neighbor exchange
            doNeighExch();
        }
        else
        {
            // ghost forward exchange
            doP2P(gh_fw_len, gh_fw_r_targets, gh_fw_s_targets, gh_fw_s_sizes, gh_fw_cyc);
        }

        // k-space pre exchange
        doP2P(k_pre_len, k_pre_r_targets, k_pre_s_targets, k_pre_s_sizes, k_pre_cyc);

        // do FFT
        doFFT();

        // k-space post exchange
        doP2P(k_post_len, k_post_r_targets, k_post_s_targets, k_post_s_sizes, k_post_cyc);

        // energy calculation
        SWM_Compute(k_energy_cyc);
        SWM_Allreduce(48, rsp_bytes, SWM_COMM_WORLD, req_vc, resp_vc, NO_BUFFER, NO_BUFFER);

        // ghost reverse exchange
        doP2P(gh_rw_len, gh_rw_r_targets, gh_rw_s_targets, gh_rw_s_sizes, gh_rw_cyc);

        // ghost fixed values exchange
        doP2P(fix_len, fix_r_targets, fix_s_targets, fix_s_sizes, fix_cyc);

        // final integration
        SWM_Compute(final_cyc);
        SWM_Allreduce(8, rsp_bytes, SWM_COMM_WORLD, req_vc, resp_vc, NO_BUFFER, NO_BUFFER);  // temperature
        SWM_Allreduce(48, rsp_bytes, SWM_COMM_WORLD, req_vc, resp_vc, NO_BUFFER, NO_BUFFER); // pressure
    }
    SWM_Finalize();
    //MM: comment assert(0);
}


#pragma GCC diagnostic push
void
LAMMPS_SWM::modelInit()
{
#include "lammps_model.h"
    double t_vol, f_vol;
    int i = 0, j = 0;

    t_vol = prd[0]/procNums[0] * prd[1]/procNums[1] * prd[2]/procNums[2];
    f_vol = pppmGrid[0]/procNums[0] * pppmGrid[1]/procNums[1] * pppmGrid[2]/procNums[2];


    get_k_params(process_id, f_vol);
    ghost_setup(GHOST_SKIN_CUTOFF, process_id, t_vol);
    k_pre_setup(FFT_SKIN_CUTOFF, process_id, f_vol);
    k_post_setup(FFT_SKIN_CUTOFF, process_id, f_vol);
    neigh_e_setup(GHOST_SKIN_CUTOFF, process_id, t_vol);

    neigh_check_average = neigh_check_avg;
    neigh_check_cyc = std::max((long)0, (long)((t_vol * ins_neigh_check_a[0] + ins_neigh_check_b[0]) * ins_neigh_check_cpi * router_freq / cpu_freq / cpu_sim_speedup + 0.5));
    for(i = 0; i < NUM_NEIGH_ALLREDUCE; i++)
    {
        neigh_end_cyc[i] = std::max((long)0, (long)((t_vol * ins_neigh_end_a[i] + ins_neigh_end_b[i]) * ins_neigh_end_cpi * router_freq / cpu_freq / cpu_sim_speedup + 0.5));
    }

    start_cyc = std::max((long)0, (long)((t_vol * ins_start_a[0] + ins_start_b[0]) * ins_start_cpi * router_freq / cpu_freq / cpu_sim_speedup + 0.5));
    k_energy_cyc = std::max((long)0, (long)((f_vol * ins_k_energy_a[0] + ins_k_energy_b[0]) * ins_k_energy_cpi * router_freq / cpu_freq / cpu_sim_speedup + 0.5));
    final_cyc = std::max((long)0, (long)((t_vol * ins_final_a[0] + ins_final_b[0]) * ins_final_cpi * router_freq / cpu_freq / cpu_sim_speedup + 0.5));


    neigh_check_count = 0;
    neigh_check_cumulative = 0.;

//     for(i = 0; i < 3; i++){
//         printf("procs: %d, %d, %d\n pppm: %lf, %lf, %lf\n", procNums[0], procNums[1], procNums[2], pppmGrid[0], pppmGrid[1], pppmGrid[2]);
//     }

//     for(j = 0; j < gh_fw_len; j++){
//         printf("%d: (%d->%d, %d->%d)\n", process_id, gh_fw_r_targets[j], process_id, process_id, gh_fw_s_targets[j]);
//     }

//     for(j = 0; j < neigh_e_len; j++){
//         printf("%d: (%d->%d, %d->%d)\n", process_id, neigh_e_r_targets[j], process_id, process_id, neigh_e_s_targets[j]);
//     }

//     if(process_id == 11){
//         for(i = 0; i < NUM_TRANSPOSE; i++){
//             printf("transpose: %d\n", i);
//             printf("len: %d, (recv proc, send_proc, send size): ", k_len[i]);
//             for(j = 0; j < k_len[i]; j++){
//                 printf("(%d, %d, %d) ", k_r_targets[i][j], k_s_targets[i][j], k_s_sizes[i][j]);
//             }
//             printf("\n\n");
//         }
//     }


//     if(process_id == 11){

//         printf("gh_fw:\n");
//         printf("len: %d, (recv proc, send_proc, send size): ", gh_fw_len);
//         for(j = 0; j < gh_fw_len; j++){
//             printf("(%d, %d, %d) ", gh_fw_r_targets[j], gh_fw_s_targets[j], gh_fw_s_sizes[j]);
//         }
//         printf("\n\n");

//         printf("gh_rw:\n");
//         printf("len: %d, (recv proc, send_proc, send size): ", gh_rw_len);
//         for(j = 0; j < gh_rw_len; j++){
//             printf("(%d, %d, %d) ", gh_rw_r_targets[j], gh_rw_s_targets[j], gh_rw_s_sizes[j]);
//         }
//         printf("\n\n");

//         printf("fix:\n");
//         printf("len: %d, (recv proc, send_proc, send size): ", fix_len);
//         for(j = 0; j < fix_len; j++){
//             printf("(%d, %d, %d) ", fix_r_targets[j], fix_s_targets[j], fix_s_sizes[j]);
//         }
//         printf("\n\n");

//         printf("neigh_e:\n");
//         printf("len: %d, (recv proc, send_proc, send size): ", neigh_e_len);
//         for(j = 0; j < neigh_e_len; j++){
//             printf("(%d, %d, %d) ", neigh_e_r_targets[j], neigh_e_s_targets[j], neigh_e_s_sizes[j]);
//         }
//         printf("\n\n");

//         printf("neigh_b:\n");
//         printf("len: %d, (recv proc, send_proc, send size): ", neigh_b_len);
//         for(j = 0; j < neigh_b_len; j++){
//             printf("(%d, %d, %d) ", neigh_b_r_targets[j], neigh_b_s_targets[j], neigh_b_s_sizes[j]);
//         }
//         printf("\n\n");

//         printf("k_pre:\n");
//         printf("len: %d, (recv proc, send_proc, send size): ", k_pre_len);
//         for(j = 0; j < k_pre_len; j++){
//             printf("(%d, %d, %d) ", k_pre_r_targets[j], k_pre_s_targets[j], k_pre_s_sizes[j]);
//         }
//         printf("\n\n");

//         printf("k_post:\n");
//         printf("len: %d, (recv proc, send_proc, send size): ", k_post_len);
//         for(j = 0; j < k_post_len; j++){
//             printf("(%d, %d, %d) ", k_post_r_targets[j], k_post_s_targets[j], k_post_s_sizes[j]);
//         }
//         printf("\n\n");

//     }
}
#pragma GCC diagnostic pop

void
LAMMPS_SWM::proc_decomposition(int n, double prd[], int procNums[])
{
    double area[3];
    double bestArea;
    double tmpArea;
    int i = 0, j = 0;

    procNums[0]= procNums[1]= procNums[2]=0;

    area[0]=prd[0]*prd[1];
    area[1]=prd[0]*prd[2];
    area[2]=prd[1]*prd[2];
    bestArea = 2*(area[0]+area[1]+area[2]);

    for(i = 1; i <= n; i++)
    {
        if(n%i == 0)
        {
            for(j = 1; j <= n/i; j++)
            {
                if(n/i%j == 0)
                {
                    tmpArea = area[0]/i/j + area[1]/i/(n/i/j) + area[2]/j/(n/i/j);
                    if(tmpArea < bestArea)
                    {
                        bestArea = tmpArea;
                        procNums[0]=i;
                        procNums[1]=j;
                        procNums[2]=n/i/j;
                    }
                }
            }
        }
    }
}
void
LAMMPS_SWM::pppm_decomposition(int n, double prd[], double pppmGrid[])
{
    double h[3];
    double err;
    int i;

    h[0] = h[1] = h[2] = 1./GEWALD;


    for(i = 0; i < 3; i++) pppmGrid[i] = int(prd[i]/h[i]) + 1;

    for(i = 0; i < 3; i++)
    {
        err = pppm_estimate_ik_error(h[i], prd[i], n, prd);
        while(err > FFT_ACCURACY)
        {
            err = pppm_estimate_ik_error(h[i], prd[i], n, prd);
            pppmGrid[i]++;
            h[i] = prd[i]/pppmGrid[i];
        }
    }


    for(i = 0; i < 3; i++) while(pppm_factorable(pppmGrid[i]) == 0) pppmGrid[i]++;

}

#pragma GCC diagnostic push
void
LAMMPS_SWM::ghost_setup(double cutoff, int rank, double t_vol)
{
#include "lammps_model.h"
    int nc[3], neigh[6];
    int i = 0, ni = 0;
    double tmp_vol, max_vol;

    gh_fw_len = 0;
    for(i = 0; i < 3; i++) nc[i] = int(cutoff / (prd[i] / procNums[i]) + 1);
    for(i = 0; i < 3; i++) gh_fw_len+=2*nc[i];

    gh_fw_r_targets = new int[gh_fw_len];
    gh_fw_s_targets = new int[gh_fw_len];
    gh_fw_s_sizes = new int[gh_fw_len];
    gh_fw_cyc = new long[gh_fw_len];

//    printf("\n Rank id %d gh_fw_len %d ", rank, gh_fw_len);
    // receive targets
    rank_to_neigh(rank, neigh);
    ni = 0;
    for(i = 0; i < nc[0]; i++)
    {
        gh_fw_r_targets[ni++] = neigh[0];
        gh_fw_r_targets[ni++] = neigh[1];
    }
    for(i = 0; i < nc[1]; i++)
    {
        gh_fw_r_targets[ni++] = neigh[2];
        gh_fw_r_targets[ni++] = neigh[3];
    }
    for(i = 0; i < nc[2]; i++)
    {
        gh_fw_r_targets[ni++] = neigh[4];
        gh_fw_r_targets[ni++] = neigh[5];
    }


    // send targets and sizes
    ni = 0;
    tmp_vol = 0;
    max_vol = (prd[1]/procNums[1]+0*cutoff)*(prd[2]/procNums[2]+0*cutoff)*cutoff;
    for(i = 0; i < nc[0]; i++)
    {
        gh_fw_s_targets[ni] = neigh[1];
        if(i < nc[0]-1)
        {
            gh_fw_s_sizes[ni] = (prd[1]/procNums[1]+0*cutoff)*(prd[2]/procNums[2]+0*cutoff)*prd[0]/procNums[0];
            tmp_vol = tmp_vol + gh_fw_s_sizes[ni];
        }
        else
        {
            gh_fw_s_sizes[ni] = max_vol - tmp_vol;
        }
        ni++;

        gh_fw_s_targets[ni] = neigh[0];
        gh_fw_s_sizes[ni] = gh_fw_s_sizes[ni-1];
        ni++;
    }

    tmp_vol = 0;
    max_vol = (prd[0]/procNums[0]+2*cutoff)*(prd[2]/procNums[2]+0*cutoff)*cutoff;
    for(i = 0; i < nc[1]; i++)
    {
        gh_fw_s_targets[ni] = neigh[3];
        if(i < nc[1]-1)
        {
            gh_fw_s_sizes[ni] = (prd[0]/procNums[0]+2*cutoff)*(prd[2]/procNums[2]+0*cutoff)*prd[1]/procNums[1];
            tmp_vol = tmp_vol + gh_fw_s_sizes[ni];
        }
        else
        {
            gh_fw_s_sizes[ni] = max_vol - tmp_vol;
        }
        ni++;

        gh_fw_s_targets[ni] = neigh[2];
        gh_fw_s_sizes[ni] = gh_fw_s_sizes[ni-1];
        ni++;
    }

    tmp_vol = 0;
    max_vol = (prd[0]/procNums[0]+2*cutoff)*(prd[1]/procNums[1]+2*cutoff)*cutoff;
    for(i = 0; i < nc[2]; i++)
    {
        gh_fw_s_targets[ni] = neigh[5];
        if(i < nc[2]-1)
        {
            gh_fw_s_sizes[ni] = (prd[0]/procNums[0]+2*cutoff)*(prd[1]/procNums[1]+2*cutoff)*prd[2]/procNums[2];
            tmp_vol = tmp_vol + gh_fw_s_sizes[ni];
        }
        else
        {
            gh_fw_s_sizes[ni] = max_vol - tmp_vol;
        }
        ni++;

        gh_fw_s_targets[ni] = neigh[4];
        gh_fw_s_sizes[ni] = gh_fw_s_sizes[ni-1];
        ni++;
    }

    // reverse
    gh_rw_len = gh_fw_len;
    gh_rw_r_targets = new int[gh_rw_len];
    gh_rw_s_targets = new int[gh_rw_len];
    gh_rw_s_sizes = new int[gh_rw_len];
    gh_rw_cyc = new long[gh_rw_len];

    ni = 0;

    for(i = gh_fw_len-2; i >= 0; i=i-2)
    {
        gh_rw_r_targets[ni] = gh_fw_r_targets[i];
        gh_rw_s_targets[ni] = gh_fw_s_targets[i];
        gh_rw_s_sizes[ni] = gh_fw_s_sizes[i];
        ni++;

        gh_rw_r_targets[ni] = gh_fw_r_targets[i+1];
        gh_rw_s_targets[ni] = gh_fw_s_targets[i+1];
        gh_rw_s_sizes[ni] = gh_fw_s_sizes[i+1];
        ni++;
    }


    // fix setup
    fix_len = gh_fw_len;
    fix_r_targets = new int[fix_len];
    fix_s_targets = new int[fix_len];
    fix_s_sizes = new int[fix_len];
    fix_cyc = new long[fix_len];

    for(i = 0; i < gh_fw_len; i++)
    {
        fix_r_targets[i] = gh_fw_r_targets[i];
        fix_s_targets[i] = gh_fw_s_targets[i];
        fix_s_sizes[i] = gh_fw_s_sizes[i];
    }


    // neigh_borders setup
    neigh_b_len = gh_fw_len;
    neigh_b_r_targets = new int[neigh_b_len];
    neigh_b_s_targets = new int[neigh_b_len];
    neigh_b_s_sizes = new int[neigh_b_len];
    neigh_b_cyc = new long[neigh_b_len];

    for(i = 0; i < gh_fw_len; i++)
    {
        neigh_b_r_targets[i] = gh_fw_r_targets[i];
        neigh_b_s_targets[i] = gh_fw_s_targets[i];
        neigh_b_s_sizes[i] = gh_fw_s_sizes[i];
    }


    // scale the sizes
    for(i = 0; i < gh_fw_len; i++)   gh_fw_s_sizes[i]   = (int)(gh_fw_s_sizes[i]   * msg_ghost_fw     + 0.5);
    for(i = 0; i < gh_rw_len; i++)   gh_rw_s_sizes[i]   = (int)(gh_rw_s_sizes[i]   * msg_ghost_rw     + 0.5);
    for(i = 0; i < fix_len; i++)     fix_s_sizes[i]     = (int)(fix_s_sizes[i]     * msg_fix          + 0.5);
    for(i = 0; i < neigh_b_len; i++) neigh_b_s_sizes[i] = (int)(neigh_b_s_sizes[i] * msg_neigh_border + 0.5);


    // instruction counts
    ni = 0;
    for(i = 0; i < nc[0]; i++)
    {
        gh_fw_cyc[ni]   = t_vol * ins_ghost_fw_a[0]        + ins_ghost_fw_b[0];
        fix_cyc[ni]     = t_vol * ins_fix_a[0]             + ins_fix_b[0];
        neigh_b_cyc[ni] = t_vol * ins_neigh_border_sr_a[0] + ins_neigh_border_sr_b[0];
        ni++;
        gh_fw_cyc[ni]   = t_vol * ins_ghost_fw_a[1]        + ins_ghost_fw_b[1];
        fix_cyc[ni]     = t_vol * ins_fix_a[1]             + ins_fix_b[1];
        neigh_b_cyc[ni] = t_vol * ins_neigh_border_sr_a[1] + ins_neigh_border_sr_b[1];
        ni++;
    }
    for(i = 0; i < nc[1]; i++)
    {
        gh_fw_cyc[ni]   = t_vol * ins_ghost_fw_a[2]        + ins_ghost_fw_b[2];
        fix_cyc[ni]     = t_vol * ins_fix_a[2]             + ins_fix_b[2];
        neigh_b_cyc[ni] = t_vol * ins_neigh_border_sr_a[2] + ins_neigh_border_sr_b[2];
        gh_fw_cyc[ni] = t_vol * ins_ghost_fw_a[2] + ins_ghost_fw_b[2];
        ni++;
        gh_fw_cyc[ni]   = t_vol * ins_ghost_fw_a[3]        + ins_ghost_fw_b[3];
        fix_cyc[ni]     = t_vol * ins_fix_a[3]             + ins_fix_b[3];
        neigh_b_cyc[ni] = t_vol * ins_neigh_border_sr_a[3] + ins_neigh_border_sr_b[3];
        ni++;
    }
    for(i = 0; i < nc[2]; i++)
    {
        gh_fw_cyc[ni]   = t_vol * ins_ghost_fw_a[4]        + ins_ghost_fw_b[4];
        fix_cyc[ni]     = t_vol * ins_fix_a[4]             + ins_fix_b[4];
        neigh_b_cyc[ni] = t_vol * ins_neigh_border_sr_a[4] + ins_neigh_border_sr_b[4];
        ni++;
        gh_fw_cyc[ni]   = t_vol * ins_ghost_fw_a[5]        + ins_ghost_fw_b[5];
        fix_cyc[ni]     = t_vol * ins_fix_a[5]             + ins_fix_b[5];
        neigh_b_cyc[ni] = t_vol * ins_neigh_border_sr_a[5] + ins_neigh_border_sr_b[5];
        ni++;
    }
    ni = 0;
    for(i = 0; i < nc[2]; i++)
    {
        gh_rw_cyc[ni++] = t_vol * ins_ghost_rw_a[4] + ins_ghost_rw_b[4];
        gh_rw_cyc[ni++] = t_vol * ins_ghost_rw_a[5] + ins_ghost_rw_b[5];
    }
    for(i = 0; i < nc[1]; i++)
    {
        gh_rw_cyc[ni++] = t_vol * ins_ghost_rw_a[2] + ins_ghost_rw_b[2];
        gh_rw_cyc[ni++] = t_vol * ins_ghost_rw_a[3] + ins_ghost_rw_b[3];
    }
    for(i = 0; i < nc[0]; i++)
    {
        gh_rw_cyc[ni++] = t_vol * ins_ghost_rw_a[0] + ins_ghost_rw_b[0];
        gh_rw_cyc[ni++] = t_vol * ins_ghost_rw_a[1] + ins_ghost_rw_b[1];
    }


    // convert instructions to cycles
    for(i = 0; i < gh_fw_len; i++)
        gh_fw_cyc[i] = std::max((long)0, (long)(gh_fw_cyc[i] * ins_ghost_fw_cpi * router_freq / cpu_freq / cpu_sim_speedup + 0.5));
    for(i = 0; i < gh_rw_len; i++)
        gh_rw_cyc[i] = std::max((long)0, (long)(gh_rw_cyc[i] * ins_ghost_rw_cpi * router_freq / cpu_freq / cpu_sim_speedup + 0.5));
    for(i = 0; i < fix_len; i++)
        fix_cyc[i] = std::max((long)0, (long)(fix_cyc[i] * ins_fix_cpi * router_freq / cpu_freq / cpu_sim_speedup + 0.5));
    for(i = 0; i < neigh_b_len; i++)
        neigh_b_cyc[i] = std::max((long)0, (long)(neigh_b_cyc[i] * ins_neigh_border_sr_cpi * router_freq / cpu_freq / cpu_sim_speedup + 0.5));

}

void
LAMMPS_SWM::k_pre_setup(double cutoff, int rank, double f_vol)
{
#include "lammps_model.h"
    int neigh[6];
    int i = 0, ni = 0;
    int hi_in, hi_out, lo_in, lo_out;
    int rs[3];
    int coord[3];

    //printf("\n proc[0] %d proc[1] %d proc[2] %d ", procNums[0], procNums[1], procNums[2]);
    //printf("\n proc[0] %f proc[1] %f proc[2] %f cutoff %f ", prd[0], prd[1], prd[2], cutoff);
    for(i = 0; i < 3; i++) 
    {
//        printf("\n cutoff %f prd[%d] %f procNums[%d] %d ", cutoff, i, prd[i], i, procNums[i]);
        assert(int(cutoff / (prd[i] / procNums[i]) + 1) == 1);
    }
    k_pre_len = 6;

    k_pre_r_targets = new int[k_pre_len];
    k_pre_s_targets = new int[k_pre_len];
    k_pre_s_sizes = new int[k_pre_len];
    k_pre_cyc = new long[k_pre_len];

    // receive targets
    rank_to_neigh(rank, neigh);
    ni = 0;
    k_pre_r_targets[ni++] = neigh[4];
    k_pre_r_targets[ni++] = neigh[5];
    k_pre_r_targets[ni++] = neigh[2];
    k_pre_r_targets[ni++] = neigh[3];
    k_pre_r_targets[ni++] = neigh[0];
    k_pre_r_targets[ni++] = neigh[1];

    // send targets
    ni = 0;
    k_pre_s_targets[ni++] = neigh[5];
    k_pre_s_targets[ni++] = neigh[4];
    k_pre_s_targets[ni++] = neigh[3];
    k_pre_s_targets[ni++] = neigh[2];
    k_pre_s_targets[ni++] = neigh[1];
    k_pre_s_targets[ni++] = neigh[0];


    // send sizes
    rank_to_xyz(rank, coord);
    ni = 0;

    // receive sizes
    hi_out = (int)(((coord[0]+1)*prd[0]/procNums[0] + cutoff/2.0) * pppmGrid[0]/prd[0] + 0.5);
    hi_in = (int)((coord[0]+1)*prd[0]/procNums[0] * pppmGrid[0]/prd[0]) - 1;
    lo_out = (int)((coord[0]*prd[0]/procNums[0] - cutoff/2.0) * pppmGrid[0]/prd[0] + 0.5);
    lo_in = (int)(coord[0]*prd[0]/procNums[0] * pppmGrid[0]/prd[0]);
    rs[0] = abs(lo_out - lo_in) + 2 + abs(hi_out - hi_in) + 2;
    hi_out = (int)(((coord[1]+1)*prd[1]/procNums[1] + cutoff/2.0) * pppmGrid[1]/prd[1] + 0.5);
    hi_in = (int)((coord[1]+1)*prd[1]/procNums[1] * pppmGrid[1]/prd[1]) - 1;
    lo_out = (int)((coord[1]*prd[1]/procNums[1] - cutoff/2.0) * pppmGrid[1]/prd[1] + 0.5);
    lo_in = (int)(coord[1]*prd[1]/procNums[1] * pppmGrid[1]/prd[1]);
    rs[1] = abs(lo_out - lo_in) + 2 + abs(hi_out - hi_in) + 2;
    hi_out = (int)(((coord[2]+1)*prd[2]/procNums[2] + cutoff/2.0) * pppmGrid[2]/prd[2] + 0.5);
    hi_in = (int)((coord[2]+1)*prd[2]/procNums[2] * pppmGrid[2]/prd[2]) - 1;
    lo_out = (int)((coord[2]*prd[2]/procNums[2] - cutoff/2.0) * pppmGrid[2]/prd[2] + 0.5);
    lo_in = (int)(coord[2]*prd[2]/procNums[2] * pppmGrid[2]/prd[2]);
    rs[2] = abs(lo_out - lo_in) + 2 + abs(hi_out - hi_in) + 2;

    // send sizes
    lo_out = (int)round(((coord[2]+1)*prd[2]/procNums[2] + cutoff/2.0) * pppmGrid[2]/prd[2] + 0.5, 10);
    lo_in = (int)round((coord[2]+1)*prd[2]/procNums[2] * pppmGrid[2]/prd[2], 10) - 1;
    hi_out = (int)round((coord[2]*prd[2]/procNums[2] - cutoff/2.0) * pppmGrid[2]/prd[2] + 0.5, 10);
    hi_in = (int)round(coord[2]*prd[2]/procNums[2] * pppmGrid[2]/prd[2], 10);
    k_pre_s_sizes[ni++] = (int)((abs(hi_out - hi_in) + 2) *
                                (((int)(pppmGrid[0]/procNums[0]*(coord[0]+1))-(int)(pppmGrid[0]/procNums[0]*coord[0])) + rs[0]) *
                                (((int)(pppmGrid[1]/procNums[1]*(coord[1]+1))-(int)(pppmGrid[1]/procNums[1]*coord[1])) + rs[1]));
    k_pre_s_sizes[ni++] = (int)((abs(lo_out - lo_in) + 2) *
                                (((int)(pppmGrid[0]/procNums[0]*(coord[0]+1))-(int)(pppmGrid[0]/procNums[0]*coord[0])) + rs[0]) *
                                (((int)(pppmGrid[1]/procNums[1]*(coord[1]+1))-(int)(pppmGrid[1]/procNums[1]*coord[1])) + rs[1]));

    assert( (int)(abs(hi_out - hi_in) + 2) <=
            ((int)pppmGrid[2]/procNums[2] + (int)(double)((int)pppmGrid[2]%procNums[2])/procNums[2]*(coord[2]+1)) );
    assert( (int)(abs(lo_out - lo_in) + 2) <=
            ((int)pppmGrid[2]/procNums[2] + (int)(double)((int)pppmGrid[2]%procNums[2])/procNums[2]*(coord[2]+1)) );

    lo_out = (int)round(((coord[1]+1)*prd[1]/procNums[1] + cutoff/2.0) * pppmGrid[1]/prd[1] + 0.5, 10);
    lo_in = (int)round((coord[1]+1)*prd[1]/procNums[1] * pppmGrid[1]/prd[1], 10) - 1;
    hi_out = (int)round((coord[1]*prd[1]/procNums[1] - cutoff/2.0) * pppmGrid[1]/prd[1] + 0.5, 10);
    hi_in = (int)round(coord[1]*prd[1]/procNums[1] * pppmGrid[1]/prd[1], 10);
    k_pre_s_sizes[ni++] = (int)((abs(hi_out - hi_in) + 2) *
                                (((int)(pppmGrid[0]/procNums[0]*(coord[0]+1))-(int)(pppmGrid[0]/procNums[0]*coord[0])) + rs[0]) *
                                ((int)(pppmGrid[2]/procNums[2]*(coord[2]+1))-(int)(pppmGrid[2]/procNums[2]*coord[2])));
    k_pre_s_sizes[ni++] = (int)((abs(lo_out - lo_in) + 2) *
                                (((int)(pppmGrid[0]/procNums[0]*(coord[0]+1))-(int)(pppmGrid[0]/procNums[0]*coord[0])) + rs[0]) *
                                ((int)(pppmGrid[2]/procNums[2]*(coord[2]+1))-(int)(pppmGrid[2]/procNums[2]*coord[2])));

    assert( (int)(abs(hi_out - hi_in) + 2) <=
            ((int)pppmGrid[1]/procNums[1] + (int)(double)((int)pppmGrid[1]%procNums[1])/procNums[1]*(coord[1]+1)) );
    assert( (int)(abs(lo_out - lo_in) + 2) <=
            ((int)pppmGrid[1]/procNums[1] + (int)(double)((int)pppmGrid[1]%procNums[1])/procNums[1]*(coord[1]+1)) );

    lo_out = (int)round(((coord[0]+1)*prd[0]/procNums[0] + cutoff/2.0) * pppmGrid[0]/prd[0] + 0.5, 10);
    lo_in = (int)round((coord[0]+1)*prd[0]/procNums[0] * pppmGrid[0]/prd[0], 10) - 1;
    hi_out = (int)round((coord[0]*prd[0]/procNums[0] - cutoff/2.0) * pppmGrid[0]/prd[0] + 0.5, 10);
    hi_in = (int)round(coord[0]*prd[0]/procNums[0] * pppmGrid[0]/prd[0], 10);
    k_pre_s_sizes[ni++] = (int)((abs(hi_out - hi_in) + 2) *
                                ((int)(pppmGrid[1]/procNums[1]*(coord[1]+1))-(int)(pppmGrid[1]/procNums[1]*coord[1])) *
                                ((int)(pppmGrid[2]/procNums[2]*(coord[2]+1))-(int)(pppmGrid[2]/procNums[2]*coord[2])));
    k_pre_s_sizes[ni++] = (int)((abs(lo_out - lo_in) + 2) *
                                ((int)(pppmGrid[1]/procNums[1]*(coord[1]+1))-(int)(pppmGrid[1]/procNums[1]*coord[1])) *
                                ((int)(pppmGrid[2]/procNums[2]*(coord[2]+1))-(int)(pppmGrid[2]/procNums[2]*coord[2])));

    assert( (int)(abs(hi_out - hi_in) + 2) <=
            ((int)pppmGrid[0]/procNums[0] + (int)(double)((int)pppmGrid[0]%procNums[0])/procNums[0]*(coord[0]+1)) );
    assert( (int)(abs(lo_out - lo_in) + 2) <=
            ((int)pppmGrid[0]/procNums[0] + (int)(double)((int)pppmGrid[0]%procNums[0])/procNums[0]*(coord[0]+1)) );

    printf("\n lo_out-lo_in %d down-val %d ", lo_out - lo_in, (int)pppmGrid[0]/procNums[0] + (int)(double)((int)pppmGrid[0]%procNums[0])/procNums[0]*(coord[0]+1));

    for(i = 0; i < k_pre_len; i++) k_pre_s_sizes[i] = int(k_pre_s_sizes[i] * msg_k_pre + 0.5);

    // cycle counts
    for(i = 0; i < k_pre_len; i++)
    {
        k_pre_cyc[i] = std::max((long)0, (long)((f_vol * ins_k_pre_a[i] + ins_k_pre_b[i]) * ins_k_pre_cpi * router_freq / cpu_freq / cpu_sim_speedup + 0.5));
    }

}

void
LAMMPS_SWM::k_post_setup(double cutoff, int rank, double f_vol)
{
#include "lammps_model.h"
    int neigh[6];
    int i, ni;
    int hi_in, hi_out, lo_in, lo_out;
    int rs[3];
    int coord[3];

    for(i = 0; i < 3; i++) assert(int(cutoff / (prd[i] / procNums[i]) + 1) == 1);
    k_post_len = 6;

    k_post_r_targets = new int[k_post_len];
    k_post_s_targets = new int[k_post_len];
    k_post_s_sizes = new int[k_post_len];
    k_post_cyc = new long[k_post_len];

    // receive targets
    rank_to_neigh(rank, neigh);
    ni = 0;
    k_post_r_targets[ni++] = neigh[0];
    k_post_r_targets[ni++] = neigh[1];
    k_post_r_targets[ni++] = neigh[2];
    k_post_r_targets[ni++] = neigh[3];
    k_post_r_targets[ni++] = neigh[4];
    k_post_r_targets[ni++] = neigh[5];

    // send targets
    ni = 0;
    k_post_s_targets[ni++] = neigh[1];
    k_post_s_targets[ni++] = neigh[0];
    k_post_s_targets[ni++] = neigh[3];
    k_post_s_targets[ni++] = neigh[2];
    k_post_s_targets[ni++] = neigh[5];
    k_post_s_targets[ni++] = neigh[4];


    // send sizes
    rank_to_xyz(rank, coord);
    ni = 0;


    // receive sizes
    hi_out = (int)(((coord[0]+1)*prd[0]/procNums[0] + cutoff/2.0) * pppmGrid[0]/prd[0] + 0.5);
    hi_in = (int)((coord[0]+1)*prd[0]/procNums[0] * pppmGrid[0]/prd[0]) - 1;
    lo_out = (int)((coord[0]*prd[0]/procNums[0] - cutoff/2.0) * pppmGrid[0]/prd[0] + 0.5);
    lo_in = (int)(coord[0]*prd[0]/procNums[0] * pppmGrid[0]/prd[0]);
    rs[0] = (abs(lo_out - lo_in) + 2 + abs(hi_out - hi_in) + 2);
    hi_out = (int)(((coord[1]+1)*prd[1]/procNums[1] + cutoff/2.0) * pppmGrid[1]/prd[1] + 0.5);
    hi_in = (int)((coord[1]+1)*prd[1]/procNums[1] * pppmGrid[1]/prd[1]) - 1;
    lo_out = (int)((coord[1]*prd[1]/procNums[1] - cutoff/2.0) * pppmGrid[1]/prd[1] + 0.5);
    lo_in = (int)(coord[1]*prd[1]/procNums[1] * pppmGrid[1]/prd[1]);
    rs[1] = (abs(lo_out - lo_in) + 2 + abs(hi_out - hi_in) + 2);
    hi_out = (int)(((coord[2]+1)*prd[2]/procNums[2] + cutoff/2.0) * pppmGrid[2]/prd[2] + 0.5);
    hi_in = (int)((coord[2]+1)*prd[2]/procNums[2] * pppmGrid[2]/prd[2]) - 1;
    lo_out = (int)((coord[2]*prd[2]/procNums[2] - cutoff/2.0) * pppmGrid[2]/prd[2] + 0.5);
    lo_in = (int)(coord[2]*prd[2]/procNums[2] * pppmGrid[2]/prd[2]);
    rs[2] = (abs(lo_out - lo_in) + 2 + abs(hi_out - hi_in) + 2);

    // send sizes
    lo_out = (int)round((((coord[0]-1)%procNums[0]+1)*prd[0]/procNums[0] + cutoff/2.0) * pppmGrid[0]/prd[0] + 0.5, 10);
    lo_in = (int)round(((coord[0]-1)%procNums[0]+1)*prd[0]/procNums[0] * pppmGrid[0]/prd[0], 10) - 1;
    hi_out = (int)round(((coord[0]+1)%procNums[0]*prd[0]/procNums[0] - cutoff/2.0) * pppmGrid[0]/prd[0] + 0.5, 10);
    hi_in = (int)round((coord[0]+1)%procNums[0]*prd[0]/procNums[0] * pppmGrid[0]/prd[0], 10);
    k_post_s_sizes[ni++] = (int)((abs(lo_out - lo_in) + 2) *
                                 ((int)(pppmGrid[1]/procNums[1]*(coord[1]+1))-(int)(pppmGrid[1]/procNums[1]*coord[1])) *
                                 ((int)(pppmGrid[2]/procNums[2]*(coord[2]+1))-(int)(pppmGrid[2]/procNums[2]*coord[2])));
    k_post_s_sizes[ni++] = (int)((abs(hi_out - hi_in) + 2) *
                                 ((int)(pppmGrid[1]/procNums[1]*(coord[1]+1))-(int)(pppmGrid[1]/procNums[1]*coord[1])) *
                                 ((int)(pppmGrid[2]/procNums[2]*(coord[2]+1))-(int)(pppmGrid[2]/procNums[2]*coord[2])));

    assert( (int)(abs(hi_out - hi_in) + 2) <=
            (pppmGrid[0]/procNums[0] + (int)(double)((int)pppmGrid[0]%procNums[0])/procNums[0]*(coord[0]+1)));
    assert( (int)(abs(lo_out - lo_in) + 2) <=
            (pppmGrid[0]/procNums[0] + (int)(double)((int)pppmGrid[0]%procNums[0])/procNums[0]*(coord[0]+1)));

    lo_out = (int)round((((coord[1]-1)%procNums[1]+1)*prd[1]/procNums[1] + cutoff/2.0) * pppmGrid[1]/prd[1] + 0.5, 10);
    lo_in = (int)round(((coord[1]-1)%procNums[1]+1)*prd[1]/procNums[1] * pppmGrid[1]/prd[1], 10) - 1;
    hi_out = (int)round(((coord[1]+1)%procNums[1]*prd[1]/procNums[1] - cutoff/2.0) * pppmGrid[1]/prd[1] + 0.5, 10);
    hi_in = (int)round((coord[1]+1)%procNums[1]*prd[1]/procNums[1] * pppmGrid[1]/prd[1], 10);
    k_post_s_sizes[ni++] = (int)((abs(lo_out - lo_in) + 2) *
                                 (((int)(pppmGrid[0]/procNums[0]*(coord[0]+1))-(int)(pppmGrid[0]/procNums[0]*coord[0])) + rs[0]) *
                                 ((int)(pppmGrid[2]/procNums[2]*(coord[2]+1))-(int)(pppmGrid[2]/procNums[2]*coord[2])));
    k_post_s_sizes[ni++] = (int)((abs(hi_out - hi_in) + 2) *
                                 (((int)(pppmGrid[0]/procNums[0]*(coord[0]+1))-(int)(pppmGrid[0]/procNums[0]*coord[0])) + rs[0]) *
                                 ((int)(pppmGrid[2]/procNums[2]*(coord[2]+1))-(int)(pppmGrid[2]/procNums[2]*coord[2])));

    assert( (int)(abs(hi_out - hi_in) + 2) <=
            (pppmGrid[1]/procNums[1] + (int)(double)((int)pppmGrid[1]%procNums[1])/procNums[1]*(coord[1]+1)));
    assert( (int)(abs(lo_out - lo_in) + 2) <=
            (pppmGrid[1]/procNums[1] + (int)(double)((int)pppmGrid[1]%procNums[1])/procNums[1]*(coord[1]+1)));

    lo_out = (int)round((((coord[2]-1)%procNums[2]+1)*prd[2]/procNums[2] + cutoff/2.0) * pppmGrid[2]/prd[2] + 0.5, 10);
    lo_in = (int)round(((coord[2]-1)%procNums[2]+1)*prd[2]/procNums[2] * pppmGrid[2]/prd[2], 10) - 1;
    hi_out = (int)round(((coord[2]+1)%procNums[2]*prd[2]/procNums[2] - cutoff/2.0) * pppmGrid[2]/prd[2] + 0.5, 10);
    hi_in = (int)round((coord[2]+1)%procNums[2]*prd[2]/procNums[2] * pppmGrid[2]/prd[2], 10);
    k_post_s_sizes[ni++] = (int)((abs(lo_out - lo_in) + 2) *
                                 (((int)(pppmGrid[0]/procNums[0]*(coord[0]+1))-(int)(pppmGrid[0]/procNums[0]*coord[0])) + rs[0]) *
                                 (((int)(pppmGrid[1]/procNums[1]*(coord[1]+1))-(int)(pppmGrid[1]/procNums[1]*coord[1])) + rs[1]));
    k_post_s_sizes[ni++] = (int)((abs(hi_out - hi_in) + 2) *
                                 (((int)(pppmGrid[0]/procNums[0]*(coord[0]+1))-(int)(pppmGrid[0]/procNums[0]*coord[0])) + rs[0]) *
                                 (((int)(pppmGrid[1]/procNums[1]*(coord[1]+1))-(int)(pppmGrid[1]/procNums[1]*coord[1])) + rs[1]));

    assert( (int)(abs(hi_out - hi_in) + 2) <=
            (pppmGrid[2]/procNums[2] + (int)(double)((int)pppmGrid[2]%procNums[2])/procNums[2]*(coord[2]+1)));
    assert( (int)(abs(lo_out - lo_in) + 2) <=
            (pppmGrid[2]/procNums[2] + (int)(double)((int)pppmGrid[2]%procNums[2])/procNums[2]*(coord[2]+1)));


    for(i = 0; i < k_post_len; i++) k_post_s_sizes[i] = int(k_post_s_sizes[i] * msg_k_post + 0.5);

    // cycle counts
    for(i = 0; i < k_post_len; i++)
    {
        k_post_cyc[i] = std::max((long)0, (long)((f_vol * ins_k_post_a[i] + ins_k_post_b[i]) * ins_k_post_cpi * router_freq / cpu_freq / cpu_sim_speedup + 0.5));
    }

}


void
LAMMPS_SWM::neigh_e_setup(double cutoff, int rank, double t_vol)
{
#include "lammps_model.h"
    int neigh[6];
    int i = 0, ni = 0;

    neigh_e_len = 0;
    rank_to_neigh(rank, neigh);

    for(i = 0; i < 6; i=i+2)
    {
        neigh_e_len++;
        if(neigh[i] != neigh[i+1]) neigh_e_len++;
    }

    neigh_e_r_targets = new int[neigh_e_len];
    neigh_e_s_targets = new int[neigh_e_len];
    neigh_e_s_sizes = new int[neigh_e_len];
    neigh_e_cyc = new long[neigh_e_len];

    // receive targets
    ni = 0;
    neigh_e_r_targets[ni++] = neigh[0];
    if(neigh[0] != neigh[1]) neigh_e_r_targets[ni++] = neigh[1];
    neigh_e_r_targets[ni++] = neigh[2];
    if(neigh[2] != neigh[3]) neigh_e_r_targets[ni++] = neigh[3];
    neigh_e_r_targets[ni++] = neigh[4];
    if(neigh[4] != neigh[5]) neigh_e_r_targets[ni++] = neigh[5];

    // send targets
    ni = 0;
    neigh_e_s_targets[ni++] = neigh[1];
    if(neigh[0] != neigh[1]) neigh_e_s_targets[ni++] = neigh[0];
    neigh_e_s_targets[ni++] = neigh[3];
    if(neigh[2] != neigh[3]) neigh_e_s_targets[ni++] = neigh[2];
    neigh_e_s_targets[ni++] = neigh[5];
    if(neigh[4] != neigh[5]) neigh_e_s_targets[ni++] = neigh[4];

    // send sizes
    ni = 0;
    neigh_e_s_sizes[ni++] = (int)(prd[1]/procNums[1])*(prd[2]/procNums[2]);
    if(neigh[0] != neigh[1]) neigh_e_s_sizes[ni++] = (int)(prd[1]/procNums[1])*(prd[2]/procNums[2]);
    neigh_e_s_sizes[ni++] = (int)(prd[0]/procNums[0])*(prd[2]/procNums[2]);
    if(neigh[2] != neigh[3]) neigh_e_s_sizes[ni++] = (int)(prd[0]/procNums[0])*(prd[2]/procNums[2]);
    neigh_e_s_sizes[ni++] = (int)(prd[0]/procNums[0])*(prd[1]/procNums[1]);
    if(neigh[4] != neigh[5]) neigh_e_s_sizes[ni++] = (int)(prd[0]/procNums[0])*(prd[1]/procNums[1]);

    for(i = 0; i < neigh_e_len; i++) neigh_e_s_sizes[i] = int(neigh_e_s_sizes[i] * msg_neigh_exch + 0.5);

    // setup cycle counts
    ni = 0;
    neigh_e_cyc[ni++] = std::max((long)0, (long)((t_vol * ins_neigh_exch_sr_a[0] + ins_neigh_exch_sr_b[0]) * ins_neigh_exch_sr_cpi * router_freq / cpu_freq / cpu_sim_speedup + 0.5));
    if(neigh[0] != neigh[1]) neigh_e_cyc[ni++] = 0;
    neigh_e_cyc[ni++] = std::max((long)0, (long)((t_vol * ins_neigh_exch_sr_a[1] + ins_neigh_exch_sr_b[1]) * ins_neigh_exch_sr_cpi * router_freq / cpu_freq / cpu_sim_speedup + 0.5));
    if(neigh[2] != neigh[3]) neigh_e_cyc[ni++] = 0;
    neigh_e_cyc[ni++] = std::max((long)0, (long)((t_vol * ins_neigh_exch_sr_a[2] + ins_neigh_exch_sr_b[2]) * ins_neigh_exch_sr_cpi * router_freq / cpu_freq / cpu_sim_speedup + 0.5));
    if(neigh[4] != neigh[5]) neigh_e_cyc[ni++] = 0;

}


double
LAMMPS_SWM::pppm_estimate_ik_error(double h, double p, int n, double prd[])
{
    double q2;
    double acons[5] = {1.0/23232.0, 7601.0/13628160.0, 143.0/69120.0, 517231.0/106536960.0, 106640677.0/11737571328.0};
    double sum = 0;
    int i;

    q2 = 19.426017 * sqrt(n*prd[0]*prd[1]*prd[2]);

    for(i = 0; i < 5; i++) sum += acons[i] * pow(h*GEWALD, 2.0*i);

    return q2*pow(h*GEWALD,5)*sqrt(GEWALD*p*sqrt(2*PI)*sum/n)/(p*p);

}


int
LAMMPS_SWM::pppm_factorable(int n)
{
    int factors[3] = {2, 3, 5};
    int i = 0, flag = 0;

    while(n > 1)
    {
        flag = 1;
        for(i = 0; i < 3; i++)
        {
            if(n%factors[i] == 0)
            {
                n = n / factors[i];
                flag = 0;
                break;
            }
        }
        if(flag == 1) return 0;
    }
    return 1;
}

void
LAMMPS_SWM::get_k_params(int rank, double f_vol)
{
#include "lammps_model.h"
    int *r_r, *s_r, *s_rs;
    int r_len, s_len;
    int *nx_in, *nx_fft, *nx_mid1, *nx_mid2;
    int i = 0, n_tr = 0, j = 0;

    r_r = new int[process_cnt];
    s_r = new int[process_cnt];
    s_rs = new int[process_cnt];
    nx_in = new int[10*process_cnt];
    nx_fft = new int[10*process_cnt];
    nx_mid1 = new int[10*process_cnt];
    nx_mid2 = new int[10*process_cnt];

    for(i = 0; i < (int)process_cnt; i++)
    {
        get_nx_in(i, &nx_in[i*10]);
        get_nx_fft(i, &nx_fft[i*10]);
        get_nx_mid1(i, &nx_mid1[i*10]);
        get_nx_mid2(i, &nx_mid2[i*10]);
    }

    n_tr = 0;
    find_overlap(nx_in, 0, nx_fft, 0, rank, r_r, &r_len, s_r, s_rs, &s_len);
    printf("\n r_len %d s_len %d rank %d ", r_len, s_len, rank);
    assert(r_len == s_len);
    assert(n_tr < NUM_TRANSPOSE);
    k_r_targets[n_tr] = new int[r_len];
    k_s_targets[n_tr] = new int[r_len];
    k_s_sizes[n_tr] = new int[r_len];
    for(i = 0; i < r_len; i++)
    {
        k_r_targets[n_tr][i] = r_r[i];
        k_s_targets[n_tr][i] = s_r[i];
        k_s_sizes[n_tr][i] = 8*s_rs[i];
    }
    k_len[n_tr] = r_len;

    n_tr++;
    find_overlap(nx_fft, 0, nx_mid1, 0, rank, r_r, &r_len, s_r, s_rs, &s_len);
    assert(r_len == s_len);
    assert(n_tr < NUM_TRANSPOSE);
    k_r_targets[n_tr] = new int[r_len];
    k_s_targets[n_tr] = new int[r_len];
    k_s_sizes[n_tr] = new int[r_len];
    for(i = 0; i < r_len; i++)
    {
        k_r_targets[n_tr][i] = r_r[i];
        k_s_targets[n_tr][i] = s_r[i];
        k_s_sizes[n_tr][i] = 16*s_rs[i];
    }
    k_len[n_tr] = r_len;

    n_tr++;
    find_overlap(nx_mid1, 2, nx_mid2, 2, rank, r_r, &r_len, s_r, s_rs, &s_len);
    assert(r_len == s_len);
    assert(n_tr < NUM_TRANSPOSE);
    k_r_targets[n_tr] = new int[r_len];
    k_s_targets[n_tr] = new int[r_len];
    k_s_sizes[n_tr] = new int[r_len];
    for(i = 0; i < r_len; i++)
    {
        k_r_targets[n_tr][i] = r_r[i];
        k_s_targets[n_tr][i] = s_r[i];
        k_s_sizes[n_tr][i] = 16*s_rs[i];
    }
    k_len[n_tr] = r_len;

    n_tr++;
    find_overlap(nx_mid2, 4, nx_fft, 4, rank, r_r, &r_len, s_r, s_rs, &s_len);
    assert(r_len == s_len);
    assert(n_tr < NUM_TRANSPOSE);
    k_r_targets[n_tr] = new int[r_len];
    k_s_targets[n_tr] = new int[r_len];
    k_s_sizes[n_tr] = new int[r_len];
    for(i = 0; i < r_len; i++)
    {
        k_r_targets[n_tr][i] = r_r[i];
        k_s_targets[n_tr][i] = s_r[i];
        k_s_sizes[n_tr][i] = 16*s_rs[i];
    }
    k_len[n_tr] = r_len;


    for(j = 0; j < 3; j++)
    {
        n_tr++;
        find_overlap(nx_fft, 0, nx_mid1, 0, rank, r_r, &r_len, s_r, s_rs, &s_len);
        assert(r_len == s_len);
        assert(n_tr < NUM_TRANSPOSE);
        k_r_targets[n_tr] = new int[r_len];
        k_s_targets[n_tr] = new int[r_len];
        k_s_sizes[n_tr] = new int[r_len];
        for(i = 0; i < r_len; i++)
        {
            k_r_targets[n_tr][i] = r_r[i];
            k_s_targets[n_tr][i] = s_r[i];
            k_s_sizes[n_tr][i] = 16*s_rs[i];
        }
        k_len[n_tr] = r_len;

        n_tr++;
        find_overlap(nx_mid1, 2, nx_mid2, 2, rank, r_r, &r_len, s_r, s_rs, &s_len);
        assert(r_len == s_len);
        assert(n_tr < NUM_TRANSPOSE);
        k_r_targets[n_tr] = new int[r_len];
        k_s_targets[n_tr] = new int[r_len];
        k_s_sizes[n_tr] = new int[r_len];
        for(i = 0; i < r_len; i++)
        {
            k_r_targets[n_tr][i] = r_r[i];
            k_s_targets[n_tr][i] = s_r[i];
            k_s_sizes[n_tr][i] = 16*s_rs[i];
        }
        k_len[n_tr] = r_len;

        n_tr++;
        find_overlap(nx_mid2, 4, nx_in, 4, rank, r_r, &r_len, s_r, s_rs, &s_len);
        assert(r_len == s_len);
        assert(n_tr < NUM_TRANSPOSE);
        k_r_targets[n_tr] = new int[r_len];
        k_s_targets[n_tr] = new int[r_len];
        k_s_sizes[n_tr] = new int[r_len];
        for(i = 0; i < r_len; i++)
        {
            k_r_targets[n_tr][i] = r_r[i];
            k_s_targets[n_tr][i] = s_r[i];
            k_s_sizes[n_tr][i] = 16*s_rs[i];
        }
        k_len[n_tr] = r_len;
    }


    delete r_r;
    delete s_r;
    delete s_rs;
    delete nx_in;
    delete nx_fft;
    delete nx_mid1;
    delete nx_mid2;

    // cycle counts
    for(i = 0; i < NUM_TRANSPOSE; i++)
    {
        k_cyc[i] = std::max((long)0, (long)((f_vol * ins_k_fft_a[i] + ins_k_fft_b[i]) * ins_k_fft_cpi * router_freq / cpu_freq / cpu_sim_speedup + 0.5));
    }
}

#pragma GCC diagnostic pop

int
LAMMPS_SWM::find_one_overlap(int a[6], int b[6], int s[3])
{
    int r[6];

    r[0] = std::max(a[0], b[0]);
    r[1] = std::min(a[1], b[1]);
    r[2] = std::max(a[2], b[2]);
    r[3] = std::min(a[3], b[3]);
    r[4] = std::max(a[4], b[4]);
    r[5] = std::min(a[5], b[5]);

    s[0] = s[1] = s[2] = 0;

    if( (r[0] > r[1]) || (r[2] > r[3]) || (r[4] > r[5]) ) return 0;

    s[0] = r[1] - r[0] + 1;
    s[1] = r[3] - r[2] + 1;
    s[2] = r[5] - r[4] + 1;

    return 1;
}

void
LAMMPS_SWM::find_overlap(int all_in[], int in_shift, int all_out[], int out_shift, int rank, int r_r[], int *r_len, int s_r[], int s_rs[], int *s_len)
{
    int i, r;
    int s[3];

    *r_len = 0;
    *s_len = 0;

    for(i = 1; i < (int)process_cnt; i++)
    {
        r = (rank + i) % process_cnt;

        if(find_one_overlap(&all_in[rank*10+in_shift], &all_out[r*10+out_shift], s))
        {
            //printf("\n s_len is %d ", *s_len);
            s_r[*s_len] = r;
            s_rs[*s_len] = s[0] * s[1] *s[2];
            (*s_len)++;
        }

        //printf("\n r is %d ", r);
        if(find_one_overlap(&all_in[r*10+in_shift], &all_out[rank*10+out_shift], s))
        {
            //printf("\n r_len is %d ", r*10+in_shift);
            r_r[*r_len] = r;
            (*r_len)++;
        }
    }
}

void
LAMMPS_SWM::get_nx_in(int rank, int nx[10])
{
    int coord[3];

    rank_to_xyz(rank, coord);

    nx[0]=int(double(coord[0]) / procNums[0] * pppmGrid[0]);
    nx[1]=int(double(coord[0]+1) / procNums[0] * pppmGrid[0]) - 1;
    nx[2]=int(double(coord[1]) / procNums[1] * pppmGrid[1]);
    nx[3]=int(double(coord[1]+1) / procNums[1] * pppmGrid[1]) - 1;
    nx[4]=int(double(coord[2]) / procNums[2] * pppmGrid[2]);
    nx[5]=int(double(coord[2]+1) / procNums[2] * pppmGrid[2]) - 1;
    nx[6]=nx[0];
    nx[7]=nx[1];
    nx[8]=nx[2];
    nx[9]=nx[3];
}

void
LAMMPS_SWM::get_nx_fft(int rank, int nx[10])
{
    int py, pz, me_y, me_z;

    if(pppmGrid[2] > process_cnt)
    {
        py = 1;
        pz = process_cnt;
    }
    else
    {
        best_2d_mapping(&py, &pz, int(pppmGrid[1]), int(pppmGrid[2]));
    }

    me_y = rank % py;
    me_z = rank / py;

    nx[0] = 0;
    nx[1] = pppmGrid[0] - 1;
    nx[2] = me_y * pppmGrid[1] / py;
    nx[3] = (me_y+1) * pppmGrid[1] / py - 1;
    nx[4] = me_z * pppmGrid[2] / pz;
    nx[5] = (me_z+1) * pppmGrid[2] / pz - 1;
    nx[6]=nx[0];
    nx[7]=nx[1];
    nx[8]=nx[2];
    nx[9]=nx[3];
}

void
LAMMPS_SWM::get_nx_mid1(int rank, int nx[10])
{
    int f1, f2;
    int ip1, ip2;

    bifactor(process_cnt, &f1, &f2);

    ip1 = rank % f1;
    ip2 = rank / f1;

    nx[0]=ip1*(int)pppmGrid[0]/f1;
    nx[1]=(ip1+1)*(int)pppmGrid[0]/f1-1;
    nx[2]=0;
    nx[3]=(int)pppmGrid[1]-1;
    nx[4]=ip2*(int)pppmGrid[2]/f2;
    nx[5]=(ip2+1)*(int)pppmGrid[2]/f2-1;
    nx[6]=nx[0];
    nx[7]=nx[1];
    nx[8]=nx[2];
    nx[9]=nx[3];
}

void
LAMMPS_SWM::get_nx_mid2(int rank, int nx[10])
{
    int f1, f2;
    int ip1, ip2;

    bifactor(process_cnt, &f1, &f2);

    ip1 = rank % f1;
    ip2 = rank / f1;

    nx[0]=ip1*(int)pppmGrid[0]/f1;
    nx[1]=(ip1+1)*(int)pppmGrid[0]/f1-1;
    nx[2]=ip2*(int)pppmGrid[1]/f2;
    nx[3]=(ip2+1)*(int)pppmGrid[1]/f2-1;
    nx[4]=0;
    nx[5]=(int)pppmGrid[2]-1;
    nx[6]=nx[0];
    nx[7]=nx[1];

    nx[8]=nx[2];
    nx[9]=nx[3];
}


void
LAMMPS_SWM::best_2d_mapping(int *px, int *py, int nx, int ny)
{
    int bestsurf, bestboxx, bestboxy;
    int boxx, boxy, surf;
    int ipx, ipy;

    bestsurf = 2 * (nx + ny);
    bestboxx = 0;
    bestboxy = 0;

    ipx = 1;

    while(ipx <= (int)process_cnt)
    {
        if(process_cnt % ipx == 0)
        {
            ipy = process_cnt / ipx;
            boxx = nx / ipx;
            if(nx % ipx) boxx++;
            boxy = ny / ipy;
            if(ny % ipy) boxy++;
            surf = boxx + boxy;
            if( (surf < bestsurf) ||
                    ( (surf == bestsurf) && (boxx*boxy > bestboxx*bestboxy) ))
            {

                bestsurf = surf;
                bestboxx = boxx;
                bestboxy = boxy;
                *px = ipx;
                *py = ipy;
            }
        }
        ipx++;
    }
}

void
LAMMPS_SWM::bifactor(int n, int *f1, int *f2)
{
    *f1 = int(sqrt(n));
    while(*f1 > 0)
    {
        *f2 = n / *f1;
        if((*f1) * (*f2) == n) return;
        (*f1)--;
    }
}

void
LAMMPS_SWM::rank_to_xyz(int rank, int coord[3])
{
    coord[0] = rank / procNums[2] / procNums[1] % procNums[0];
    coord[1] = rank / procNums[2] % procNums[1];
    coord[2] = rank % procNums[2];
}

int
LAMMPS_SWM::xyz_to_rank(int coord[3])
{
    int mods[3];
    int i = 0;
    for(i = 0; i < 3; i++)
    {
        mods[i] = coord[i]%procNums[i];
        while(mods[i] < 0) mods[i] += procNums[i];
    }
    return mods[0]*procNums[1]*procNums[2] + mods[1]*procNums[2] + mods[2];
}

void
LAMMPS_SWM::rank_to_neigh(int rank, int neighs[6])
{
    int coord[3];
    int tmp_coord[3];
    rank_to_xyz(rank, coord);

    tmp_coord[0]=coord[0]+1;
    tmp_coord[1]=coord[1];
    tmp_coord[2]=coord[2];
    neighs[0] = xyz_to_rank(tmp_coord);
    tmp_coord[0]=coord[0]-1;
    tmp_coord[1]=coord[1];
    tmp_coord[2]=coord[2];
    neighs[1] = xyz_to_rank(tmp_coord);
    tmp_coord[0]=coord[0];
    tmp_coord[1]=coord[1]+1;
    tmp_coord[2]=coord[2];
    neighs[2] = xyz_to_rank(tmp_coord);
    tmp_coord[0]=coord[0];
    tmp_coord[1]=coord[1]-1;
    tmp_coord[2]=coord[2];
    neighs[3] = xyz_to_rank(tmp_coord);
    tmp_coord[0]=coord[0];
    tmp_coord[1]=coord[1];
    tmp_coord[2]=coord[2]+1;
    neighs[4] = xyz_to_rank(tmp_coord);
    tmp_coord[0]=coord[0];
    tmp_coord[1]=coord[1];
    tmp_coord[2]=coord[2]-1;
    neighs[5] = xyz_to_rank(tmp_coord);
}


//DLL_POSTAMBLE(LAMMPS_SWM)
