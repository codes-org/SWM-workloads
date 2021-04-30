#ifndef _LAMMPS_
#define _LAMMPS_

#include <list>
#include <boost/property_tree/ptree.hpp>
#include <iostream>	// For printing message sizes
#include "swm-include.h"
// Internal LAMMPS paramenters
// Skin cutoff for ghost neighbor exchange (on comm)
#define GHOST_SKIN_CUTOFF 12.0
// Skin cutoff for fft neighbor exchange (on commgrid)
#define FFT_SKIN_CUTOFF 2.0
// Number of atoms in a basic block
#define N_ATOMS_BASE 32000
// Neighbor check after NEIGH_DELAY, then every NEIGH_EVERY
#define NEIGH_DELAY 5
#define NEIGH_EVERY 1
// Dimensions of the basic block
#define XLO_BASE (-27.5)
#define XHI_BASE (27.5)
#define YLO_BASE (-38.5)
#define YHI_BASE (38.5)
#define ZLO_BASE (-36.3646)
#define ZHI_BASE (36.3615)
// lammps factors for determining required decomposition
#define GEWALD 0.243177
#define FFT_ACCURACY 0.033206
// number of transposes in fft
#define NUM_TRANSPOSE 13
// number of allreduces at the end of neighbor exchange
#define NUM_NEIGH_ALLREDUCE 5


#define PI 3.14159265358979323846


class LAMMPS_SWM 
{
public:

    LAMMPS_SWM(

        //SWMUserIF* user_if,
        boost::property_tree::ptree cfg,
        void**& generic_ptrs
    );
    ~LAMMPS_SWM();
    void call();

protected:
    uint32_t x_rep;      // number of replicas in X dimension
    uint32_t y_rep;      // number of replicas in Y dimension
    uint32_t z_rep;      // number of replicas in Z dimension
    uint32_t num_timesteps; // number of time steps to simulate
    uint32_t req_vc;        // request vc
    uint32_t resp_vc;       // response vc
    double router_freq;     // router frequency in Hz
    double cpu_freq;        // CPU frequency in Hz
    double cpu_sim_speedup; // simulation speedup factor (makes CPU faster)
    uint32_t rsp_bytes;
    uint32_t process_cnt; // MM addition
    int process_id; //MM addition

private:
    double prd[3];
    double pppmGrid[3];
    int procNums[3];

    int *k_r_targets[NUM_TRANSPOSE];
    int *k_s_targets[NUM_TRANSPOSE];
    int *k_s_sizes[NUM_TRANSPOSE];
    long k_cyc[NUM_TRANSPOSE];
    int k_len[NUM_TRANSPOSE];

    int *gh_fw_r_targets;
    int *gh_fw_s_targets;
    int *gh_fw_s_sizes;
    long *gh_fw_cyc;
    int gh_fw_len;

    int *gh_rw_r_targets;
    int *gh_rw_s_targets;
    int *gh_rw_s_sizes;
    long *gh_rw_cyc;
    int gh_rw_len;

    int *k_pre_r_targets;
    int *k_pre_s_targets;
    int *k_pre_s_sizes;
    long *k_pre_cyc;
    int k_pre_len;

    int *k_post_r_targets;
    int *k_post_s_targets;
    int *k_post_s_sizes;
    long *k_post_cyc;
    int k_post_len;

    int *fix_r_targets;
    int *fix_s_targets;
    int *fix_s_sizes;
    long *fix_cyc;
    int fix_len;

    int *neigh_e_r_targets;
    int *neigh_e_s_targets;
    int *neigh_e_s_sizes;
    long *neigh_e_cyc;
    int neigh_e_len;

    int *neigh_b_r_targets;
    int *neigh_b_s_targets;
    int *neigh_b_s_sizes;
    long *neigh_b_cyc;
    int neigh_b_len;

    long neigh_check_cyc;
    double neigh_check_average;
    double neigh_check_cumulative;
    int neigh_check_count;
    long neigh_end_cyc[NUM_NEIGH_ALLREDUCE];

    long start_cyc;
    long k_energy_cyc;
    long final_cyc;

    void modelInit();
    void doP2P(int len, int *r_targets, int *s_targets, int *s_sizes, long *cyc_cnt);
    void doNeighExch();
    void doFFT();
    bool neigh_check();

    // process decomposition
    void proc_decomposition(int n, double prd[], int procNums[]);

    // PPPM decomposition
    void pppm_decomposition(int n, double prd[], double pppmGrid[]);
    double pppm_estimate_ik_error(double h, double prd, int n, double all_prd[]);
    int pppm_factorable(int n);

    // neighbor comm setup
    void ghost_setup(double cutoff, int rank, double t_vol);
    void k_pre_setup(double cutoff, int rank, double f_vol);
    void k_post_setup(double cutoff, int rank, double f_vol);
    void neigh_e_setup(double cutoff, int rank, double t_vol);

    // k space paramenters
    void get_k_params(int rank, double f_vol);
    void get_nx_in(int rank, int nx[10]);
    void get_nx_fft(int rank, int nx[10]);
    void get_nx_mid1(int rank, int nx[10]);
    void get_nx_mid2(int rank, int nx[10]);

    int find_one_overlap(int a[6], int b[6], int s[3]);
    void find_overlap(int all_in[], int in_shift, int all_out[], int out_shift, int rank, int r_r[], int *r_len, int s_r[], int s_rs[], int *s_len);

    void best_2d_mapping(int *px, int *py, int nx, int ny);
    void bifactor(int n, int *f1, int *f2);

    void rank_to_xyz(int rank, int coord[3]);
    int xyz_to_rank(int coord[3]);
    void rank_to_neigh(int rank, int neighs[6]);
};


#endif
