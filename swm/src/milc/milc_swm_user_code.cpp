#include "milc_swm_user_code.h"
#include <climits>

template <typename T>
std::vector<T> boost_ptree_array_to_std_vector(boost::property_tree::ptree const& pt, boost::property_tree::ptree::key_type const& key, std::vector<T> def, bool disallow_empty_arrays=true)
{
  std::vector<T> r;
  for (auto& item : pt.get_child(key)) {
    r.push_back(item.second.get_value<T>());
  }
  if(r.empty()) {
    std::cerr << "\nERROR: boost_ptree_array_to_std_vector matched a tag (" << key << ") that is not a valid array!\n" << std::endl;
    assert(r.empty() == false);
  }
  return r;
}

MilcSWMUserCode::MilcSWMUserCode(
    boost::property_tree::ptree cfg,
    void**& generic_ptrs
) :
    process_cnt(cfg.get<uint32_t>("jobs.size", 1)),
    iteration_cnt(cfg.get<uint32_t>("jobs.cfg.iteration_cnt", 1)),
    compute_delay(cfg.get<uint32_t>("jobs.cfg.compute_delay", 1)),
    dimension_cnt(cfg.get<uint32_t>("jobs.cfg.dimension_cnt",0)),
    msg_size(cfg.get<uint32_t>("jobs.cfg.msg_size", 0)),
    dimension_sizes(boost_ptree_array_to_std_vector<uint32_t>(cfg,"jobs.cfg.dimension_sizes", {0})),
    max_dimension_distance(cfg.get<uint32_t>("jobs.cfg.max_dimension_distance",0)),
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
}

void
MilcSWMUserCode::xlat_pid_to_coords(uint32_t pid,
				    std::vector<uint32_t>& coords)
{
  coords.clear();

  uint32_t dim_div = 1;
  for(uint32_t dim_idx=0; dim_idx<dimension_cnt; dim_idx++) {
    uint32_t pid_coord_in_dim = (pid / dim_div) % dimension_sizes[dim_idx];
    dim_div *= dimension_sizes[dim_idx];
    coords.push_back(pid_coord_in_dim);
  }
}

void
MilcSWMUserCode::xlat_coords_to_pid(std::vector<uint32_t> coords,
				    uint32_t& pid)
{
  pid=0;
  uint32_t dim_mult = 1;
  for(uint32_t dim_idx=0; dim_idx<dimension_cnt; dim_idx++) {
    pid += coords[dim_idx] * dim_mult;
    dim_mult *= dimension_sizes[dim_idx];
  }
}

std::string
MilcSWMUserCode::get_neighbor_string(
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
MilcSWMUserCode::derive_neighbors_recurse(
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
        //std::cout << "neighbor_string is " << neighbor_string << ", regexd_string is " << regexed_string << std::endl;

        neighbors.push_back( std::make_tuple(neighbor_pid, neighbor_string) );

        return;
    }
    else if(dimension_to_vary == dimension_cnt)
    {
        if(accumulated_dimension_distance > 0)
        {
            uint32_t neighbor_pid;
            xlat_coords_to_pid(coords, neighbor_pid);

            std::string neighbor_string = get_neighbor_string(process_id, neighbor_pid);
            //std::cout << "neighbor_string is " << neighbor_string << ", regexd_string is " << regexed_string << std::endl;

            neighbors.push_back( std::make_tuple(neighbor_pid, neighbor_string) );
        }
        return;
    }

    //negative
    coords_copy = coords;
    if(coords_copy[dimension_to_vary] == 0) {
      coords_copy[dimension_to_vary] = (dimension_sizes[dimension_to_vary] -1);
    } else {
      coords_copy[dimension_to_vary] = (coords_copy[dimension_to_vary] -1);
    }
    derive_neighbors_recurse(coords_copy, neighbors,
			     dimension_to_vary+1,
			     accumulated_dimension_distance+1);

    //positive
    coords_copy = coords;
    if(coords_copy[dimension_to_vary]==(dimension_sizes[dimension_to_vary]-1)) {
      coords_copy[dimension_to_vary] = 0;
    } else {
      coords_copy[dimension_to_vary] = (coords_copy[dimension_to_vary] +1);
    }
    derive_neighbors_recurse(coords_copy, neighbors,
			     dimension_to_vary+1,
			     accumulated_dimension_distance+1);

    //none
    coords_copy = coords;
    derive_neighbors_recurse(coords_copy, neighbors,
			     dimension_to_vary+1,
			     accumulated_dimension_distance);

}

void
MilcSWMUserCode::call()
{
  /*
  if(process_id == 0) {   //lets print every pid in coords and back again
    std::vector<uint32_t> coords;
    uint32_t pid_again;
    for(uint32_t pid=0; pid<process_cnt; pid++) {
      coords.clear();
      pid_again=0;
      xlat_pid_to_coords(pid, coords);
    }
  }
  */
  struct swm_app_data app_data = {
    .final_iteration = std::min(static_cast<int>(iteration_cnt - 1), INT_MAX),
  };
  SWM_Pass_app_data(&app_data);

  std::vector<uint32_t> my_coords;
  //std::vector<uint32_t> neighbor_pids;
  xlat_pid_to_coords(process_id, my_coords);

  derive_neighbors_recurse(my_coords, neighbors);

  /*
    if(process_id == 0) {
    std::cout << "neighbors of pid " << process_id << " are: ";
    for(size_t neighbors_idx=0; neighbors_idx<neighbors.size(); neighbors_idx++) {
    std::cout << " " << std::get<0>(neighbors[neighbors_idx]) << "," << std::get<1>(neighbors[neighbors_idx]);
    }
    std::cout << "\n";
    }
  */

  uint32_t* send_handles = new uint32_t[neighbors.size()];
  uint32_t* recv_handles = new uint32_t[neighbors.size()];
  uint32_t neighbors_size=neighbors.size();
  uint32_t pkt_rsp_bytes = 0;

  for(uint32_t iter=0; iter<iteration_cnt; iter++) {
    if (process_id == 0) {
      printf("MILC: Iteration %d/%d\n",iter+1,iteration_cnt);
    }

    //shuffle the neighbors
    if(randomize_communication_order) {
      std::random_shuffle(neighbors.begin(), neighbors.end());
    }

    // Dslash on even sites, then odd sites
    for(int eo=0; eo<2; eo++) {
      // start recv's
      for(size_t neighbor_idx=0;neighbor_idx<neighbors.size();neighbor_idx++) {
	SWM_Irecv(
		  std::get<0>(neighbors[neighbor_idx]),
		  SWM_COMM_WORLD,
		  std::get<0>(neighbors[neighbor_idx]),
		  NO_BUFFER,
		  &(recv_handles[neighbor_idx])
		  );
      }
      // send to each neighbor
      for(size_t neighbor_idx=0;neighbor_idx<neighbors.size();neighbor_idx++) {
	SWM_Isend(std::get<0>(neighbors[neighbor_idx]),
		  SWM_COMM_WORLD,
		  process_id,
		  request_vc,
		  response_vc,
		  NO_BUFFER,
		  msg_size,
		  pkt_rsp_bytes,
		  &(send_handles[neighbor_idx]),
		  req_rt,
		  rsp_rt);
      }
      // local work
      SWM_Compute(compute_delay);
      SWM_Waitall(neighbors.size(), send_handles);
      SWM_Waitall(neighbors.size(), recv_handles);
    }

    for(int ar=0; ar<2; ar++) {
      SWM_Allreduce(sizeof(double), // payload
		    pkt_rsp_bytes, // pkt_rsp_bytes
		    SWM_COMM_WORLD,
		    request_vc,
		    response_vc,
		    NO_BUFFER,
		    NO_BUFFER,
		    AUTO,
		    NULL,
		    req_rt,
		    rsp_rt
		    );
    }

    SWM_Mark_Iteration(iter);

  }
  SWM_Finalize();
}
