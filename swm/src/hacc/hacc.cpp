
#include "hacc_swm_user_code.h"

HACCSWMUserCode::HACCSWMUserCode(
        SWMUserIF* user_if,
        boost::property_tree::ptree cfg,
        void**& generic_ptrs
        ) :
    SWMUserCode(user_if),
    user_if(user_if),
    request_vc(cfg.get<uint32_t>("request_vc", 0)), 
    response_vc(cfg.get<uint32_t>("response_vc", 4)),
    pkt_rsp_bytes(cfg.get<uint32_t>("pkt_rsp_bytes", 0)),
    gen_cfg_filename(cfg.get<std::string>("gen_cfg_filename")),
    cfg(cfg)
{
    //cerr << "Got filename: " << gen_cfg_filename << endl;

    //here we should parse the json pointed to by gen_cfg_filename and do whatever with it...
    std::ifstream gen_cfg_file;

    // in case there are environment paths in the variable...
    gen_cfg_filename = expand_environment_variables(gen_cfg_filename);

    gen_cfg_file.open(gen_cfg_filename);
    ASSERT(gen_cfg_file.is_open(), "Could not open gen_cfg_file: " << gen_cfg_filename << std::endl);

    //boost::property_tree::ptree gen_cfg;
    std::stringstream ss;
    std::string line;

    while(getline(gen_cfg_file, line)) {
        ss << line;
    }

    boost::property_tree::read_json(ss, gen_cfg);

    ng = gen_cfg.get<int>("ng");
    nranks = gen_cfg.get<int>("nranks"); //8;

    assert(sscanf(gen_cfg.get<std::string>("rank_shape_3d").c_str(), "(%d, %d, %d)", &(rank_shape_3d[0]), &(rank_shape_3d[1]), &(rank_shape_3d[2])) == 3);
    assert(sscanf(gen_cfg.get<std::string>("rank_shape_2d_x").c_str(), "(%d, %d, %d)", &(rank_shape_2d_x[0]), &(rank_shape_2d_x[1]), &(rank_shape_2d_x[2])) == 3);
    assert(sscanf(gen_cfg.get<std::string>("rank_shape_2d_y").c_str(), "(%d, %d, %d)", &(rank_shape_2d_y[0]), &(rank_shape_2d_y[1]), &(rank_shape_2d_y[2])) == 3);
    assert(sscanf(gen_cfg.get<std::string>("rank_shape_2d_z").c_str(), "(%d, %d, %d)", &(rank_shape_2d_z[0]), &(rank_shape_2d_z[1]), &(rank_shape_2d_z[2])) == 3);

    box_length = gen_cfg.get<double>("box_length"); //96.1458;

    /*
    printf("ng: %d\n", ng);
    printf("nranks: %d\n", nranks);
    printf("box_length: %g\n", box_length);

    for(int i=0; i<3; i++) printf("rank_shape_3d[%d]: %d\n", i, rank_shape_3d[i]);
    for(int i=0; i<3; i++) printf("rank_shape_2d_x[%d]: %d\n", i, rank_shape_2d_x[i]);
    for(int i=0; i<3; i++) printf("rank_shape_2d_y[%d]: %d\n", i, rank_shape_2d_y[i]);
    for(int i=0; i<3; i++) printf("rank_shape_2d_z[%d]: %d\n", i, rank_shape_2d_z[i]);
    */
    
    gen_cfg_file.close();
}

void
HACCSWMUserCode::call() {

    /*
    const int ng = 8;
    const int nranks = 8;
    const double box_length = 96.1458;
    const int rank_shape_3d  [3] = {2, 2, 2};
    const int rank_shape_2d_x[3] = {1, 4, 2};
    const int rank_shape_2d_y[3] = {4, 1, 2};
    const int rank_shape_2d_z[3] = {4, 2, 1};
    */

    //double box_length = 19613.75;


    // Perf model parameters
    
    const double ninteractions_per_rank_mean  = 1e10;
    const double ninteractions_per_rank_delta = 0.01;
    const double ninteractions_per_rank_per_wallsecond = 1e9;
    
    const double buffer_copy_MBps = 1000.0;
    const double fft_work_per_second = 1e9;

    bool enable_hacc_fft = cfg.get<bool>("enable_hacc_fft",true);
    bool enable_hacc_exchange = cfg.get<bool>("enable_hacc_exchange",true);
    bool enable_hacc_checksum = cfg.get<bool>("enable_hacc_checksum",true);

    // Configuration for this run
    HaccConfig config(
            ng,
            box_length,
            process_cnt,
            process_id,
            rank_shape_3d,
            rank_shape_2d_x,
            rank_shape_2d_y,
            rank_shape_2d_z,
            request_vc,
            response_vc,
            pkt_rsp_bytes
            );

    // Assemble timestep model
    timestep = new HaccTimestep (
            user_if,
            &done_to_child,
            config,

            ninteractions_per_rank_mean,
            ninteractions_per_rank_delta,
            ninteractions_per_rank_per_wallsecond,

            buffer_copy_MBps,
            fft_work_per_second,

            enable_hacc_fft,
            enable_hacc_exchange,
            enable_hacc_checksum
            );

    // ====================================================================
    // Code below will go to call() method of SWMUserCode subclass
    // ====================================================================

    // Go!
    if (enable_contexts)
        while(1) {
            (*timestep)(); //timestep.do_steps();
            if(done_to_child) break;
            else yield();
        }

    else
        (*timestep)(); //timestep.do_steps();

    SWM_Finalize();
//    assert(0);
    
}

DLL_POSTAMBLE(HACCSWMUserCode)
