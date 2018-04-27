/*
 * =====================================================================================
 *
 *       Filename:  nearest_neighbor_swm_user_code.h
 *
 *    Description:
 *
 *        Version:  1.0
 *        Created:  07/10/2013 11:11:02 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  John Thompson (FPD), john.d.thompson@intel.com
 *        Company:  Intel
 *
 * =====================================================================================
 */

#ifndef _NEAREST_NEIGHBOR_TEMPLATE_USER_CODE_
#define _NEAREST_NEIGHBOR_TEMPLATE_USER_CODE_

#include <boost/property_tree/ptree.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

#include <string>
#include <iostream>
#include <random>
#include <algorithm>
#include <vector>
#include <regex>

#include "swm-include.h"

using namespace std;

typedef std::tuple<uint32_t, std::string> neighbor_tuple;

typedef uint32_t RoutingType;

/*
struct neighbor {
    std::string neighbor_name;  //this will be of the form ([+,-]([0-9]+))+, for example +x-y+z is +0-1+2 and +x-z is +0-2
    uint32_t neighbor_pid;
};
*/

// MM: This msg traffic code is part of the base swm class

typedef struct msg_traffic_desc {

  //SWMMessageAttribute attribute;

  RoutingType msg_req_routing_type;
  RoutingType msg_rsp_routing_type;
  RoutingType pkt_rsp_routing_type;

  std::vector<uint32_t> msg_req_bytess;
  std::vector<uint32_t> msg_rsp_bytess;
  std::vector<uint32_t> pkt_rsp_bytess;

  uint32_t msg_req_bytes;
  uint32_t msg_rsp_bytes;
  uint32_t pkt_rsp_bytes;


  std::vector<uint32_t> msg_req_vcs;
  std::vector<uint32_t> msg_rsp_vcs;
  std::vector<uint32_t> pkt_rsp_vcs;

  uint32_t msg_req_vc;
  uint32_t msg_rsp_vc;
  uint32_t pkt_rsp_vc;
  
  /*
  #ifdef FABSIM_EMULATION
  stl_l2_encoding l2_encoding;
  std::vector<uint32_t> dlid_xors;
  
  uint32_t dlid_xor;
  #endif
  */
  
} msg_traffic_desc;


struct msg_traffic_set
{
  std::vector < std::tuple < uint32_t, msg_traffic_desc > > set;
  uint32_t msg_parts_sum;
  std::string name;
  std::string regex_string;

public:
  msg_traffic_set(std::string name, std::string regex_string)
  :
  msg_parts_sum(0),
    name(name),
    regex_string(regex_string)
  {}
};

std::vector<msg_traffic_set*> msg_traffic_def_vector; // MM addition: normally part of base class

class NearestNeighborSWMUserCode
{

public:

    NearestNeighborSWMUserCode(
        boost::property_tree::ptree cfg,
        void**& generic_ptrs
    );

    void xlat_pid_to_coords(uint32_t pid, std::vector<uint32_t>& coords);
    void xlat_coords_to_pid(std::vector<uint32_t> coords, uint32_t& pid);
    std::string get_neighbor_string(uint32_t my_pid, uint32_t neighbor_pid);

    void derive_neighbors_recurse(
        std::vector<uint32_t> coords,
        //std::vector<uint32_t>& neighbors,
        std::vector<neighbor_tuple>& neighbors,
        uint32_t dimension_to_vary=0,
        uint32_t accumulated_hamming_distance=0
    );

    void call();

protected:
    uint32_t process_cnt; // MM addition
    uint32_t request_vc;
    uint32_t response_vc;
    uint32_t message_size;
    uint32_t iteration_cnt; //MM addition
    uint32_t noop_cnt; //MM addition
    uint32_t compute_delay; //MM addition
    uint32_t msg_size; // MM addition
    int process_id; //MM addition
    int req_rt; // MM addition
    int rsp_rt; // MM addition
    
    uint32_t dimension_cnt;
    std::vector<uint32_t> dimension_sizes;
    uint32_t max_dimension_distance;
    bool synchronous;
    uint32_t iterations_per_sync;
    bool randomize_communication_order;

    //RoutingType req_rt;
    //RoutingType rsp_rt;

    //std::vector<uint32_t> neighbors; //we'll encode each dimension as 0, 1, or 2 in 2-bit aligned quantities
    std::vector<neighbor_tuple> neighbors; //we'll encode each dimension as 0, 1, or 2 in 2-bit aligned quantities

};

#endif
