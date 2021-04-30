/*
 * =====================================================================================
 *
 *       Filename:  nekbone_swm_user_code.h
 *
 *    Description:
 *
 *        Version:  1.0
 *        Created:  2013 Spet 11
 *       Revision:  none
 *       Compiler:  what works
 *
 *         Author:  dvillenx
 *        Company:  Intel
 *
 * =====================================================================================
 */

#ifndef _NEKBONE_TEMPLATE_USER_CODE_
#define _NEKBONE_TEMPLATE_USER_CODE_

#include <stdio.h>
#include <boost/property_tree/ptree.hpp>
#include <stdlib.h>
#include <iostream>	// For debugging output
#include "swm-include.h"


#include "triplet.h"
//#include "sizedata.h"


typedef int Err_t; //Set to zero for success.

void debug_PrintSize(int in_rank=0, bool in_outputOnlyForRank0=false);

class NEKBONESWMUserCode
{
public:
    ~NEKBONESWMUserCode();

    Err_t create();
    Err_t destroy();
    Err_t run();

    NEKBONESWMUserCode(
        boost::property_tree::ptree cfg,
        void**& generic_ptrs
    );
    void call();

protected:
    unsigned int request_vc;
    unsigned int response_vc;
    uint32_t pkt_rsp_bytes;
    uint32_t rcube;
    uint32_t ecube;
    int req_rt;
    int rsp_rt;

    /*MM additions */
    uint32_t process_cnt;
    int process_id;
    double cpu_freq;        // CPU frequency in Hz
    double cpu_sim_speedup; // simulation speedup factor (makes CPU faster)
private:
    //struct neighborLoad does the following
    //  For the current rank this->mpiRank, it haa to send "load" bytes
    //  to the neighboring rank indicated bny "rid".
    //Each rank can have at most (27-1) neighbors.
    struct NeighborLoad
    {
        Idz rid; //ID of the neighboring rank
        unsigned long load; //The number of bytes to send
    };



    int sz_nloads;  //The number of entry in nloads.
    //NeighborLoad nloads[SD::NeighborCount]; //In a cubic lattice each rank can have at most 26 neighbors.
    NeighborLoad nloads[26]; //In a cubic lattice each rank can have at most 26 neighbors.
    unsigned int largest_load_bytes; //The largest byte cound of a load in this->nloads().

    //===== Methods
    Err_t Initial_checks();
    Err_t makeMesh(int in_polyOrder);
    Err_t make_neighbors_loads();

    Err_t clearMesh();

    Err_t conjugateGradient();
    Err_t nek_gsop(const char * in_text); //Note that in_text is not use.  It is just for documentation.
    Err_t nek_glsc3();

    //===== Disabled Methods
    //NEKBONESWMUserCode( const NEKBONESWMUserCode & rhs);
    //NEKBONESWMUserCode & operator=(const NEKBONESWMUserCode & rhs);
    //
    //This specify the disribution of MPI ranks along the 3 main axis.
    //For example, (Rx,Ry,Rz) puts down
    //      Rx ranks along the X-axis
    //      Ry ranks along the Y-axis
    //      Rz ranks along the Z-axis
    // The total number of ranks in use will be Rtotal= Rx*Ry*Rz
    uint32_t Rx;
    uint32_t Ry;
    uint32_t Rz;
    uint32_t Rtotal; //Total number of MPI rank in all

    //This specify the distribution of elements within each and all ranks.
    //For example, (Ex,Ey,Ez) puts down
    //      Ex ranks along the X-axis
    //      Ey ranks along the Y-axis
    //      Ez ranks along the Z-axis
    // The total number of element in a single rank is Etotal= Ex*Ey*Ez
    uint32_t Ex;
    uint32_t Ey;
    uint32_t Ez;
    uint32_t Etotal;  //Total number of elements per rank

    //This gives the overall total element count in use
    uint32_t GlobalElementCount;

    //This specify the polynomial sequence to use.
    //For example, for (Pbegin, Pend, Pstep)=(9,13,3), then
    //  the sequence generated will be the same as
    //      for(int p=Pbegin; p!=Pend; p+=Pstep){...}
    //  where p is the polynomial order in use
    //  So p will 9,then 12.int
    //IMPORTANT Note:
    // For a given polynomial order P, NEKbone assumes that a given
    // element will be partitioned by putting P+1 degree-of-freedom (DOF)
    // along each main axis of the element.
    // So each and all elements end up with a total of (P+1)^3 DOFs each,
    // uniformaly distributed in the element.
    // For the workload specification, that ends up
    // setting (Pbegin, Pend, Pstep) = (8, 12, 3) which corresponds
    // to NEKbone::(nx0,nxN,nxD)=(9,12,3)
    uint32_t Pbegin;
    uint32_t Pend;
    uint32_t Pstep;


    //This specifies the maximum number of CG iteration to be done
    uint32_t CGcount;

    //Some global constants
    uint32_t NeighborCount; //In a cubic lattice, each cube can have at most 26 neighbors.
    uint32_t ByteSizeOf1DOF; // 8 bytes for 1 DOF --> sizeof(double)

    int check_finalize; //Only do MPI_Finalize if zero.
    /*MM changes: making this compatible with process_id type */
    int mpiRank; //The rank of the current process
    //It is also the rank ID number for the current process inside this->globPartR.
    //PARTITIONS:
    //
    //The local partitions tell how the components are split within one set:
    //      partR   --> rank distribution within the main object
    //      partE   --> element distribution within one rank
    //      partP   --> DOF distribution within one element
    //
    //The global partitions tell how a component are split over the entire object
    //      globpartR --> How all the ranks split the entire object
    //      globpartE --> How all the elements split the entire object
    //      globpartP --> How all the DOFs split the entire object
    //
    //  so   globpartR < globpartE < globpartP

    Triplet partR; //The partition along the 3 main axis for the ranks
    Triplet partE; //The partition along the 3 main axis for the element within a single rank.
    Triplet partP; //The partition along the 3 main axis for the DOF within a single element

    Triplet globPartR; //Partition over the entire mesh for only the ranks, i.e. = this->partR.
    Triplet globPartE; //Partition over the entire mesh for all the elements
    Triplet globPartP; //Partition over the entire mesh for all the DOFS.
};


#endif //_NEKBONE_TEMPLATE_USER_CODE_
