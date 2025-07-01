
#include "nearest_neighbor_swm_user_code.h"
#include "boost_ptree_array_to_std_vector.h"
extern uint64_t global_cycle;

std::vector<msg_traffic_set*> msg_traffic_def_vector; // MM addition: normally part of base class

std::string GetFirstMatch(std::string lookup_name)
{

    for(size_t s=0; s<msg_traffic_def_vector.size(); s++) {
        std::regex re(msg_traffic_def_vector[s]->regex_string);
        if(std::regex_match(lookup_name,re))
        {
        return msg_traffic_def_vector[s]->name;
        }
    }
    return "";
    //assert(0);
}

NearestNeighborSWMUserCode::NearestNeighborSWMUserCode(
    boost::property_tree::ptree cfg,
    void**& generic_ptrs
) :
    process_cnt(cfg.get<uint32_t>("jobs.size", 1)),
    iteration_cnt(cfg.get<uint32_t>("jobs.cfg.iteration_cnt", 1)),
    noop_cnt(cfg.get<uint32_t>("jobs.cfg.noop_cnt", 1)),
    compute_delay(cfg.get<uint32_t>("jobs.cfg.compute_delay", 1)),
    dimension_cnt(cfg.get<uint32_t>("jobs.cfg.dimension_cnt",0)),
    msg_size(cfg.get<uint32_t>("jobs.cfg.msg_size", 0)),
    dimension_sizes(boost_ptree_array_to_std_vector<uint32_t>(cfg,"jobs.cfg.dimension_sizes", {0})),
                max_dimension_distance(cfg.get<uint32_t>("jobs.cfg.max_dimension_distance",0)),
                synchronous(cfg.get<bool>("jobs.cfg.synchronous",true)),
                iterations_per_sync(cfg.get<uint32_t>("jobs.cfg.iterations_per_sync",1)),
                randomize_communication_order(cfg.get<bool>("jobs.cfg.randomize_communication_order",false))
{
    process_id = *((int*)generic_ptrs[0]);
    assert(dimension_sizes.size() == dimension_cnt);

    size_t dim_product = 1;
    for(size_t dim_i = 0; dim_i < dimension_sizes.size(); dim_i++)
    {
        dim_product *= dimension_sizes[dim_i];
    }
    //std::cout << "dim_product is " << dim_product << " and process_cnt is " << process_cnt << std::endl;
    assert(dim_product == process_cnt);

    req_rt = AUTOMATIC;
    rsp_rt = AUTOMATIC;

    if (synchronous == false)
        printf("SWM Nearest Neighbor - Warning: configuring 'synchronous == false' currently generates no traffic\n");

    if (randomize_communication_order && process_id == 0)
    {
        printf("SWM Nearest Neighbor - Warning: configuring 'randomize_communication_order == true' will generate nondeterminstic results.\n");
    }
}

void
NearestNeighborSWMUserCode::xlat_pid_to_coords(
    uint32_t pid,
    std::vector<uint32_t>& coords
)
{

    coords.clear();

    uint32_t dim_div = 1;
    for(uint32_t dim_idx=0; dim_idx<dimension_cnt; dim_idx++)
    {
        uint32_t pid_coord_in_dim = (pid / dim_div) % dimension_sizes[dim_idx];
        dim_div *= dimension_sizes[dim_idx];
        coords.push_back(pid_coord_in_dim);
    }
}

void
NearestNeighborSWMUserCode::xlat_coords_to_pid(
    std::vector<uint32_t> coords,
    uint32_t& pid
)
{

    pid=0;


    /*std::cout << "xlat_coords_to_pid on coords ";
    for(size_t coords_idx=0; coords_idx<coords.size(); coords_idx++) {
    std::cout << " " << coords[coords_idx];
    }
    std::cout << endl;
    */

    uint32_t dim_mult = 1;
    for(uint32_t dim_idx=0; dim_idx<dimension_cnt; dim_idx++)
    {
        pid += coords[dim_idx] * dim_mult;
        dim_mult *= dimension_sizes[dim_idx];
    }
}

std::string
NearestNeighborSWMUserCode::get_neighbor_string(
    uint32_t my_pid,
    uint32_t neighbor_pid
)
{
    std::vector<uint32_t> my_coords;
    std::vector<uint32_t> neighbor_coords;

    xlat_pid_to_coords(my_pid, my_coords);
    xlat_pid_to_coords(neighbor_pid, neighbor_coords);

    assert(my_coords.size() == neighbor_coords.size());

    std::ostringstream oss;

    for(size_t c=0; c<my_coords.size(); c++)
    {

        if(my_coords[c] != neighbor_coords[c])
        {
            if(my_coords[c] == 0)
            {
                if(neighbor_coords[c] == (my_coords[c]+1))
                {
                    oss << "p" << c;
                }
                else if(neighbor_coords[c] == (dimension_sizes[c]-1))
                {
                    oss << "m" << c;
                }
                else
                {
                    assert(0);
                }
            }
            else
            {
                if(neighbor_coords[c] == ((my_coords[c] + 1) % dimension_sizes[c]))
                {
                    oss << "p" << c;
                }
                else if(neighbor_coords[c] == (my_coords[c] - 1))
                {
                    oss << "m" << c;
                }
                else
                {
                    assert(0);
                }
            }
        }
    }

    return oss.str();
}

void
NearestNeighborSWMUserCode::derive_neighbors_recurse(
    std::vector<uint32_t> coords,
    std::vector<neighbor_tuple>& neighbors,
    uint32_t dimension_to_vary,
    uint32_t accumulated_dimension_distance
)
{
    std::vector<uint32_t> coords_copy;
    //uint32_t accumulated_dimension_distance_copy;
    coords_copy.resize(coords.size());

    if(accumulated_dimension_distance == max_dimension_distance)
    {
        uint32_t neighbor_pid;
        xlat_coords_to_pid(coords, neighbor_pid);

        std::string neighbor_string = get_neighbor_string(process_id, neighbor_pid);
        std::string regexed_string = GetFirstMatch(neighbor_string);
        //std::cout << "neighbor_string is " << neighbor_string << ", regexd_string is " << regexed_string << std::endl;

        neighbors.push_back( std::make_tuple(neighbor_pid,regexed_string) );

        return;
    }
    else if(dimension_to_vary == dimension_cnt)
    {
        if(accumulated_dimension_distance > 0)
        {
            uint32_t neighbor_pid;
            xlat_coords_to_pid(coords, neighbor_pid);

            std::string neighbor_string = get_neighbor_string(process_id, neighbor_pid);
            std::string regexed_string = GetFirstMatch(neighbor_string);
            //std::cout << "neighbor_string is " << neighbor_string << ", regexd_string is " << regexed_string << std::endl;

            neighbors.push_back( std::make_tuple(neighbor_pid,regexed_string) );
        }
        return;
    }

    //negative
    coords_copy = coords;
    if(coords_copy[dimension_to_vary] == 0)
    {
        coords_copy[dimension_to_vary] = (dimension_sizes[dimension_to_vary] -1);
    }
    else
    {
        coords_copy[dimension_to_vary] = (coords_copy[dimension_to_vary] -1);
    }

    derive_neighbors_recurse(
        coords_copy,
        neighbors,
        dimension_to_vary+1,
        accumulated_dimension_distance+1
    );


    //none
    coords_copy = coords;

    derive_neighbors_recurse(
        coords_copy,
        neighbors,
        dimension_to_vary+1,
        accumulated_dimension_distance
    );


    //positive
    coords_copy = coords;
    if(coords_copy[dimension_to_vary] == (dimension_sizes[dimension_to_vary] -1))
    {
        coords_copy[dimension_to_vary] = 0;
    }
    else
    {
        coords_copy[dimension_to_vary] = (coords_copy[dimension_to_vary] +1);
    }

    derive_neighbors_recurse(
        coords_copy,
        neighbors,
        dimension_to_vary+1,
        accumulated_dimension_distance+1
    );

}

void
NearestNeighborSWMUserCode::call()
{


    if(process_id == 0) {   //lets print every pid in coords and back again
    std::vector<uint32_t> coords;
        uint32_t pid_again;
        for(uint32_t pid=0; pid<process_cnt; pid++) {
            coords.clear();
            pid_again=0;
            xlat_pid_to_coords(pid, coords);
//            std::cout << "pid " << pid << " has coords.size " << coords.size() << " ";
//            for(size_t i=0; i<coords.size(); i++) {
//            std::cout << " " << coords[i];
//            }
//            std::cout << "; which have pid ";
//            xlat_coords_to_pid(coords, pid_again);
//            std::cout << pid_again << endl;
        }
    }


    std::vector<uint32_t> my_coords;
    std::vector<uint32_t> neighbor_pids;
    xlat_pid_to_coords(process_id, my_coords);

    derive_neighbors_recurse(my_coords, neighbors);

    /*
    if(process_id == 0)
    {
        std::cout << "neighbors of pid " << process_id << " are: ";
        for(size_t neighbors_idx=0; neighbors_idx<neighbors.size(); neighbors_idx++) {
            std::cout << " " << std::get<0>(neighbors[neighbors_idx]) << "," << std::get<1>(neighbors[neighbors_idx]);
        }
        std::cout << "\n";
    }
    */

    // uint32_t* send_handles = NULL;
    // uint32_t* recv_handles = NULL;
    uint32_t* all_handles = NULL;
    size_t num_handles_per_sync = 0;

    if(synchronous)
    {
        // send_handles = new uint32_t[neighbors.size()*iterations_per_sync];
        // recv_handles = new uint32_t[neighbors.size()*iterations_per_sync];
        all_handles = new uint32_t[neighbors.size()*iterations_per_sync * 2];
        num_handles_per_sync = neighbors.size() *  iterations_per_sync * 2;
    }

    uint32_t iter_before_sync = 0;
    uint32_t neighbors_size=neighbors.size();
    uint32_t pkt_rsp_bytes = 0;

    for(uint32_t iter=0; iter<iteration_cnt; iter++)
    {
        if(process_id == 0)
            printf("Nearest Neighbor Starting Iteration %d / %d\n",iter,iteration_cnt);

        //shuffle the neighbors
        /*NM: This is a source of nondeterminism depending on how std::rand() is implemented
              even between two sequential runs. Maybe a possible source of unmatched isends/irecv in optimistic mode*/
        if(randomize_communication_order)
        {
            std::random_shuffle(neighbors.begin(), neighbors.end());
        }

        //send to each neighbor
        for(size_t neighbor_idx=0; neighbor_idx<neighbors.size(); neighbor_idx++)
        {
            if (synchronous) {
                //send/recv pair that we'll later wait on
                SWM_Irecv(
                    std::get<0>(neighbors[neighbor_idx]),
                    SWM_COMM_WORLD,
                    std::get<0>(neighbors[neighbor_idx]),
                    NO_BUFFER,
                    &(all_handles[neighbor_idx+iter_before_sync*neighbors_size+(neighbors_size*iterations_per_sync*0)])
                );

                SWM_Isend(
                    std::get<0>(neighbors[neighbor_idx]),
                    SWM_COMM_WORLD,
                    process_id,
                    0, // MM additions
                    1, // MM additions
                    NO_BUFFER,
                    msg_size, //msg_desc.msg_req_bytes,
                    pkt_rsp_bytes, //msg_desc.pkt_rsp_bytes,
                    &(all_handles[neighbor_idx+iter_before_sync*neighbors_size+(neighbors_size*iterations_per_sync*1)]),
                    0,
                    0
                );
                for(uint32_t noop=0; noop<noop_cnt; noop++)
                {
                    assert(0); // Does CODES have an equivalent of a NOOP?
                    //SWM_Noop();
                }
            }
            else
            {
                //fire and forget
	            // purposes of the paper lets just use Isend/Irecv
	            assert(0);
	            /*
                SWM_Synthetic(
                    std::get<0>(neighbors[neighbor_idx]),  //dst
                    0,
                    0,
                    0,
                    msg_size, //msg_desc.msg_req_bytes,
                    0, //msg_desc.msg_rsp_bytes,
                    pkt_rsp_bytes, //msg_desc.pkt_rsp_bytes,
                    req_rt, //msg_desc.msg_req_routing_type,
                    rsp_rt, //msg_desc.msg_rsp_routing_type,
                    rsp_rt, //msg_desc.pkt_rsp_routing_type,
                    NULL,
                    NULL, //msg_desc.attribute
#ifdef FABSIM_EMULATION
                    , msg_desc.l2_encoding
#endif
                );
                for(uint32_t noop=0; noop<noop_cnt; noop++)
                {
                    SWM_Noop();
                }
		        */
            }
        }
        if (synchronous) {
            iter_before_sync++;
            if(iter_before_sync == iterations_per_sync || iter == iteration_cnt-1 )
            {
                //std::cout << "begin wait at time: " << global_cycle << std::endl;
                SWM_Waitall(
                    num_handles_per_sync,
                    all_handles
                );
                iter_before_sync = 0;
                //std::cout << "end wait at time: " << global_cycle << std::endl;
            }
        }
        else
        {
            if (compute_delay)
            {
                SWM_Compute(compute_delay);
            }
        }
    }
    SWM_Finalize();
}

