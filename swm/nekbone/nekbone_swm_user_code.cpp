#include "nekbone_swm_user_code.h"

#include <vector>
#include <stdio.h>
#include <string.h> //memset

#include "cubiclattice.h"


#define CHECKERR if(err) break;

#define MPI_SUCCESS (1);
typedef BUF_TYPE SWM_Request;

NEKBONESWMUserCode::NEKBONESWMUserCode(
    SWMUserIF* user_if,
    boost::property_tree::ptree cfg,
    void**& generic_ptrs
)
    :
    AppBaseSWMUserCode(user_if, cfg, "nekbone"),
    check_finalize(0), mpiRank(0),
    partR(), partE(), partP(), globPartR(), globPartE(), globPartP()
{

    msg_traffic_desc msg_desc;
    GetMsgDetails(&msg_desc);

    request_vc  = msg_desc.msg_req_vc;
    response_vc = msg_desc.msg_rsp_vc;
    pkt_rsp_bytes = msg_desc.pkt_rsp_bytes;


    assert(fabsim_cfg::get_only<uint32_t>(cfg, {"Rx"}, Rx));
    assert(fabsim_cfg::get_only<uint32_t>(cfg, {"Ry"}, Ry));
    assert(fabsim_cfg::get_only<uint32_t>(cfg, {"Rz"}, Rz));

    assert(fabsim_cfg::get_only<uint32_t>(cfg, {"Ex"}, Ex));
    assert(fabsim_cfg::get_only<uint32_t>(cfg, {"Ey"}, Ey));
    assert(fabsim_cfg::get_only<uint32_t>(cfg, {"Rx"}, Ez));

    assert(fabsim_cfg::get_only<uint32_t>(cfg, {"Pbegin"}, Pbegin));
    assert(fabsim_cfg::get_only<uint32_t>(cfg, {"Pend"}, Pend));
    assert(fabsim_cfg::get_only<uint32_t>(cfg, {"Pstep"}, Pstep));
    assert(fabsim_cfg::get_only<uint32_t>(cfg, {"CGcount"}, CGcount));
    assert(fabsim_cfg::get_only<uint32_t>(cfg, {"NeighborCount"}, NeighborCount));
    assert(fabsim_cfg::get_only<uint32_t>(cfg, {"ByteSizeOf1DOF"}, ByteSizeOf1DOF));

    RoutingTypeEnumParser routing_type_parser;
    boost::optional<std::string> req_rt_optional = cfg.get_optional<std::string>("req_rt");
    boost::optional<std::string> rsp_rt_optional = cfg.get_optional<std::string>("rsp_rt");

    if(req_rt_optional)
    {
        req_rt = routing_type_parser.ParseSomeEnum(req_rt_optional.get());
    }
    else
    {
        req_rt = AUTOMATIC;
    }

    if(rsp_rt_optional)
    {
        rsp_rt = routing_type_parser.ParseSomeEnum(rsp_rt_optional.get());
    }
    else
    {
        rsp_rt = AUTOMATIC;
    }




    //assert((rcube*rcube*rcube) == (ecube*ecube*ecube) && (rcube*rcube*rcube) == process_cnt);

    //Total number of ranks
    Rtotal = Rx*Ry*Rz;
    std::cout << "RX: " << Rx << " | RY: " << Ry << " | RZ: " << Rz << " total=RX*RY*RZ: " << Rtotal << " vs. process_cnt: " << process_cnt << std::endl;

    assert(Rtotal == process_cnt);

    //Number of elements per rank (nelt)
    Etotal = Ex*Ey*Ez;

    //Total problem size = total ranks * elements/rank
    GlobalElementCount = Rtotal*Etotal;

    Err_t err = create();
    assert(err == 0);
}


NEKBONESWMUserCode::~NEKBONESWMUserCode()
{
    Err_t err=destroy();
    assert(err == 0);
}

Err_t NEKBONESWMUserCode::create()
{
    Err_t err=0;

    //triplet.h::Idz assumes an 8-byte long.
    //If that is not the case, one run the risk of an overflow
    //when calculating IDs in cubiclattice.h
    //So flag this as an error for now.
    assert(sizeof(long) == 8);

    err = Initial_checks();

    assert(err == 0);
    return err;
}

Err_t NEKBONESWMUserCode::destroy()
{
    Err_t err=0;
    /*
       if(check_finalize>0){
    //MPI_Finalize();
    check_finalize=0;
    }
    */

    mpiRank = 0;

    err = clearMesh();

    assert(err == 0);

    //Do not delete/reset these
    //process_id;
    //request_vc;
    //response_vc

    return err;
}

void NEKBONESWMUserCode::call()
{
    Err_t err = run();
    assert (err == 0);
    SWM_Finalize();
}


Err_t NEKBONESWMUserCode::run()
{
    Err_t err=0;

    for(unsigned polyO=Pbegin; polyO<Pend; polyO+=Pstep)
    {
        if(mpiRank==0) printf("%d> polyO=%d\n", __LINE__, polyO);

        //NEKbone loop over element/rank removed-->for(G->nelt = G->iel0; G->nelt <= G->ielN; G->nelt += G->ielD){
        //Use sizedata.h::E(x|y|z) to change the element distribution within 1 rank.

        err = makeMesh(polyO);
        assert(err == 0);

        err = makeMesh(polyO);
        assert(err == 0);

        err = nek_gsop("on c");
        assert(err == 0);

        err = nek_gsop("on f");
        assert(err == 0);

        err = conjugateGradient();
        assert(err == 0);

        //mpierr = MPI_Barrier(MPI_COMM_WORLD); CHECKMPIERR;
        SWM_Barrier(
            SWM_COMM_WORLD,
            request_vc,
            response_vc,
            NO_BUFFER,
            AUTO,
            NULL,
            req_rt,
            rsp_rt
        );

        err = conjugateGradient();
        assert(err == 0);

        err = clearMesh();
        assert(err == 0);

        //NEKbone loop over element/rank removed-->}

    }

    assert(err == 0);

    return err;
}

Err_t  NEKBONESWMUserCode::Initial_checks()
{
    Err_t err=0;

    assert(!(Rx <=0 || Ry<=0 || Rz<=0));
    assert(!(Ex <=0 || Ey<=0 || Ez<=0));
    assert(!(Pbegin >= Pend || Pstep<=0));

    unsigned numtasks=0;
    //mpierr = MPI_Comm_size(MPI_COMM_WORLD, &numtasks); CHECKMPIERR;
    numtasks = process_cnt;
    assert(numtasks == Rtotal);
    //if(numtasks != Rtotal) err=__LINE__; CHECKERR;

    //mpierr = MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank); CHECKMPIERR;
    mpiRank = process_id;

    assert(!( Rtotal <= mpiRank ));

    //debug_PrintSize(mpiRank, true);

    return err;
}

Err_t NEKBONESWMUserCode::makeMesh(int in_polyOrder)
{
    Err_t err=0;
    partR = Triplet(Rx, Ry, Rz);
    partE = Triplet(Ex, Ey, Ez);
    partP = Triplet(in_polyOrder, in_polyOrder, in_polyOrder);

    globPartR = partR;
    globPartE = partE;
    globPartP = partP;

    globPartE *= globPartR;
    globPartP *= globPartE;

    err = make_neighbors_loads();

    assert(err == 0);

    return err;
}

Err_t NEKBONESWMUserCode::make_neighbors_loads()
{

#   define DEBUG_PRINT
    Err_t err=0;
    Triplet rankplet = CUBICLAT::id2tripletZ(globPartR, mpiRank);

    sz_nloads=0;
    largest_load_bytes=0;
    int t=-1; //offset in nloads
    Triplet neighbor;

    //=====The first 6 direction --> face-to-face
    neighbor = rankplet;
    neighbor += Triplet(1,0,0);
    if(neighbor.isinLattice(globPartR))
    {
        unsigned long u = partE.c * partP.c + 1;
        unsigned long v = partE.b * partP.b + 1;

        unsigned long dof = u * v * ByteSizeOf1DOF;
        if(largest_load_bytes<dof) largest_load_bytes=dof;
        Idz d = CUBICLAT::tripletZ2id(globPartR, neighbor);
        ++t;
        nloads[t].rid = d;
        nloads[t].load = dof;
        DEBUG_PRINT;
    }

    neighbor = rankplet;
    neighbor += Triplet(-1,0,0);
    if(neighbor.isinLattice(globPartR))
    {
        unsigned long u = partE.c * partP.c + 1;
        unsigned long v = partE.b * partP.b + 1;

        unsigned long dof = u * v * ByteSizeOf1DOF;
        if(largest_load_bytes<dof) largest_load_bytes=dof;
        Idz d = CUBICLAT::tripletZ2id(globPartR, neighbor);
        ++t;
        nloads[t].rid = d;
        nloads[t].load = dof;
        DEBUG_PRINT;
    }

    neighbor = rankplet;
    neighbor += Triplet(0,1,0);
    if(neighbor.isinLattice(globPartR))
    {
        unsigned long u = partE.c * partP.c + 1;
        unsigned long v = partE.a * partP.a + 1;

        unsigned long dof = u * v * ByteSizeOf1DOF;
        if(largest_load_bytes<dof) largest_load_bytes=dof;
        Idz d = CUBICLAT::tripletZ2id(globPartR, neighbor);
        ++t;
        nloads[t].rid = d;
        nloads[t].load = dof;
        DEBUG_PRINT;
    }

    neighbor = rankplet;
    neighbor += Triplet(0,-1,0);
    if(neighbor.isinLattice(globPartR))
    {
        unsigned long u = partE.c * partP.c + 1;
        unsigned long v = partE.a * partP.a + 1;

        unsigned long dof = u * v * ByteSizeOf1DOF;
        if(largest_load_bytes<dof) largest_load_bytes=dof;
        Idz d = CUBICLAT::tripletZ2id(globPartR, neighbor);
        ++t;
        nloads[t].rid = d;
        nloads[t].load = dof;
        DEBUG_PRINT;
    }

    neighbor = rankplet;
    neighbor += Triplet(0,0,1);
    if(neighbor.isinLattice(globPartR))
    {
        unsigned long u = partE.a * partP.a + 1;
        unsigned long v = partE.b * partP.b + 1;

        unsigned long dof = u * v * ByteSizeOf1DOF;
        if(largest_load_bytes<dof) largest_load_bytes=dof;
        Idz d = CUBICLAT::tripletZ2id(globPartR, neighbor);
        ++t;
        nloads[t].rid = d;
        nloads[t].load = dof;
        DEBUG_PRINT;
    }

    neighbor = rankplet;
    neighbor += Triplet(0,0,-1);
    if(neighbor.isinLattice(globPartR))
    {
        unsigned long u = partE.a * partP.a + 1;
        unsigned long v = partE.b * partP.b + 1;

        unsigned long dof = u * v * ByteSizeOf1DOF;
        if(largest_load_bytes<dof) largest_load_bytes=dof;
        Idz d = CUBICLAT::tripletZ2id(globPartR, neighbor);
        ++t;
        nloads[t].rid = d;
        nloads[t].load = dof;
        DEBUG_PRINT;
    }

    //===== The 8 directions connected corner-to-corner --> 1 DOF shared
    neighbor = rankplet;
    neighbor += Triplet(1,1,-1);
    if(neighbor.isinLattice(globPartR))
    {
        unsigned long dof = 1 * ByteSizeOf1DOF;
        if(largest_load_bytes<dof) largest_load_bytes=dof;
        Idz d = CUBICLAT::tripletZ2id(globPartR, neighbor);
        ++t;
        nloads[t].rid = d;
        nloads[t].load = dof;
        DEBUG_PRINT;
    }

    neighbor = rankplet;
    neighbor += Triplet(1,-1,-1);
    if(neighbor.isinLattice(globPartR))
    {
        unsigned long dof = 1 * ByteSizeOf1DOF;
        if(largest_load_bytes<dof) largest_load_bytes=dof;
        Idz d = CUBICLAT::tripletZ2id(globPartR, neighbor);
        ++t;
        nloads[t].rid = d;
        nloads[t].load = dof;
        DEBUG_PRINT;
    }

    neighbor = rankplet;
    neighbor += Triplet(-1,-1,-1);
    if(neighbor.isinLattice(globPartR))
    {
        unsigned long dof = 1 * ByteSizeOf1DOF;
        if(largest_load_bytes<dof) largest_load_bytes=dof;
        Idz d = CUBICLAT::tripletZ2id(globPartR, neighbor);
        ++t;
        nloads[t].rid = d;
        nloads[t].load = dof;
        DEBUG_PRINT;
    }

    neighbor = rankplet;
    neighbor += Triplet(-1,1,-1);
    if(neighbor.isinLattice(globPartR))
    {
        unsigned long dof = 1 * ByteSizeOf1DOF;
        if(largest_load_bytes<dof) largest_load_bytes=dof;
        Idz d = CUBICLAT::tripletZ2id(globPartR, neighbor);
        ++t;
        nloads[t].rid = d;
        nloads[t].load = dof;
        DEBUG_PRINT;
    }

    neighbor = rankplet;
    neighbor += Triplet(1,1,1);
    if(neighbor.isinLattice(globPartR))
    {
        unsigned long dof = 1 * ByteSizeOf1DOF;
        if(largest_load_bytes<dof) largest_load_bytes=dof;
        Idz d = CUBICLAT::tripletZ2id(globPartR, neighbor);
        ++t;
        nloads[t].rid = d;
        nloads[t].load = dof;
        DEBUG_PRINT;
    }

    neighbor = rankplet;
    neighbor += Triplet(1,-1,1);
    if(neighbor.isinLattice(globPartR))
    {
        unsigned long dof = 1 * ByteSizeOf1DOF;
        if(largest_load_bytes<dof) largest_load_bytes=dof;
        Idz d = CUBICLAT::tripletZ2id(globPartR, neighbor);
        ++t;
        nloads[t].rid = d;
        nloads[t].load = dof;
        DEBUG_PRINT;
    }

    neighbor = rankplet;
    neighbor += Triplet(-1,-1,1);
    if(neighbor.isinLattice(globPartR))
    {
        unsigned long dof = 1 * ByteSizeOf1DOF;
        if(largest_load_bytes<dof) largest_load_bytes=dof;
        Idz d = CUBICLAT::tripletZ2id(globPartR, neighbor);
        ++t;
        nloads[t].rid = d;
        nloads[t].load = dof;
        DEBUG_PRINT;
    }

    neighbor = rankplet;
    neighbor += Triplet(-1,1,1);
    if(neighbor.isinLattice(globPartR))
    {
        unsigned long dof = 1 * ByteSizeOf1DOF;
        if(largest_load_bytes<dof) largest_load_bytes=dof;
        Idz d = CUBICLAT::tripletZ2id(globPartR, neighbor);
        ++t;
        nloads[t].rid = d;
        nloads[t].load = dof;
        DEBUG_PRINT;
    }

    //===== The 12 neighbors on a edge
    neighbor = rankplet;
    neighbor += Triplet(0,-1,-1);
    if(neighbor.isinLattice(globPartR))
    {
        unsigned long u = partE.a * partP.a + 1;

        unsigned long dof = u * ByteSizeOf1DOF;
        if(largest_load_bytes<dof) largest_load_bytes=dof;
        Idz d = CUBICLAT::tripletZ2id(globPartR, neighbor);
        ++t;
        nloads[t].rid = d;
        nloads[t].load = dof;
        DEBUG_PRINT;
    }
    neighbor = rankplet;
    neighbor += Triplet(0,1,-1);
    if(neighbor.isinLattice(globPartR))
    {
        unsigned long u = partE.a * partP.a + 1;

        unsigned long dof = u * ByteSizeOf1DOF;
        if(largest_load_bytes<dof) largest_load_bytes=dof;
        Idz d = CUBICLAT::tripletZ2id(globPartR, neighbor);
        ++t;
        nloads[t].rid = d;
        nloads[t].load = dof;
        DEBUG_PRINT;
    }
    neighbor = rankplet;
    neighbor += Triplet(0,-1,1);
    if(neighbor.isinLattice(globPartR))
    {
        unsigned long u = partE.a * partP.a + 1;

        unsigned long dof = u * ByteSizeOf1DOF;
        if(largest_load_bytes<dof) largest_load_bytes=dof;
        Idz d = CUBICLAT::tripletZ2id(globPartR, neighbor);
        ++t;
        nloads[t].rid = d;
        nloads[t].load = dof;
        DEBUG_PRINT;
    }
    neighbor = rankplet;
    neighbor += Triplet(0,1,1);
    if(neighbor.isinLattice(globPartR))
    {
        unsigned long u = partE.a * partP.a + 1;

        unsigned long dof = u * ByteSizeOf1DOF;
        if(largest_load_bytes<dof) largest_load_bytes=dof;
        Idz d = CUBICLAT::tripletZ2id(globPartR, neighbor);
        ++t;
        nloads[t].rid = d;
        nloads[t].load = dof;
        DEBUG_PRINT;
    }

    neighbor = rankplet;
    neighbor += Triplet(1,0,-1);
    if(neighbor.isinLattice(globPartR))
    {
        unsigned long u = partE.b * partP.b + 1;

        unsigned long dof = u * ByteSizeOf1DOF;
        if(largest_load_bytes<dof) largest_load_bytes=dof;
        Idz d = CUBICLAT::tripletZ2id(globPartR, neighbor);
        ++t;
        nloads[t].rid = d;
        nloads[t].load = dof;
        DEBUG_PRINT;
    }
    neighbor = rankplet;
    neighbor += Triplet(1,0,1);
    if(neighbor.isinLattice(globPartR))
    {
        unsigned long u = partE.b * partP.b + 1;

        unsigned long dof = u * ByteSizeOf1DOF;
        if(largest_load_bytes<dof) largest_load_bytes=dof;
        Idz d = CUBICLAT::tripletZ2id(globPartR, neighbor);
        ++t;
        nloads[t].rid = d;
        nloads[t].load = dof;
        DEBUG_PRINT;
    }
    neighbor = rankplet;
    neighbor += Triplet(-1,0,1);
    if(neighbor.isinLattice(globPartR))
    {
        unsigned long u = partE.b * partP.b + 1;

        unsigned long dof = u * ByteSizeOf1DOF;
        if(largest_load_bytes<dof) largest_load_bytes=dof;
        Idz d = CUBICLAT::tripletZ2id(globPartR, neighbor);
        ++t;
        nloads[t].rid = d;
        nloads[t].load = dof;
        DEBUG_PRINT;
    }
    neighbor = rankplet;
    neighbor += Triplet(-1,0,-1);
    if(neighbor.isinLattice(globPartR))
    {
        unsigned long u = partE.b * partP.b + 1;

        unsigned long dof = u * ByteSizeOf1DOF;
        if(largest_load_bytes<dof) largest_load_bytes=dof;
        Idz d = CUBICLAT::tripletZ2id(globPartR, neighbor);
        ++t;
        nloads[t].rid = d;
        nloads[t].load = dof;
        DEBUG_PRINT;
    }

    neighbor = rankplet;
    neighbor += Triplet(-1,-1,0);
    if(neighbor.isinLattice(globPartR))
    {
        unsigned long u = partE.c * partP.c + 1;

        unsigned long dof = u * ByteSizeOf1DOF;
        if(largest_load_bytes<dof) largest_load_bytes=dof;
        Idz d = CUBICLAT::tripletZ2id(globPartR, neighbor);
        ++t;
        nloads[t].rid = d;
        nloads[t].load = dof;
        DEBUG_PRINT;
    }
    neighbor = rankplet;
    neighbor += Triplet(-1,1,0);
    if(neighbor.isinLattice(globPartR))
    {
        unsigned long u = partE.c * partP.c + 1;

        unsigned long dof = u * ByteSizeOf1DOF;
        if(largest_load_bytes<dof) largest_load_bytes=dof;
        Idz d = CUBICLAT::tripletZ2id(globPartR, neighbor);
        ++t;
        nloads[t].rid = d;
        nloads[t].load = dof;
        DEBUG_PRINT;
    }
    neighbor = rankplet;
    neighbor += Triplet(1,-1,0);
    if(neighbor.isinLattice(globPartR))
    {
        unsigned long u = partE.c * partP.c + 1;

        unsigned long dof = u * ByteSizeOf1DOF;
        if(largest_load_bytes<dof) largest_load_bytes=dof;
        Idz d = CUBICLAT::tripletZ2id(globPartR, neighbor);
        ++t;
        nloads[t].rid = d;
        nloads[t].load = dof;
        DEBUG_PRINT;
    }
    neighbor = rankplet;
    neighbor += Triplet(1,1,0);
    if(neighbor.isinLattice(globPartR))
    {
        unsigned long u = partE.c * partP.c + 1;

        unsigned long dof = u * ByteSizeOf1DOF;
        if(largest_load_bytes<dof) largest_load_bytes=dof;
        Idz d = CUBICLAT::tripletZ2id(globPartR, neighbor);
        ++t;
        nloads[t].rid = d;
        nloads[t].load = dof;
        DEBUG_PRINT;
    }

    //===== Finishing up
    sz_nloads = t+1;

    return err;
}

Err_t NEKBONESWMUserCode::clearMesh()
{
    Err_t err=0;
    while(!err)
    {
        memset(    (void*)&partR,0,sizeof(Triplet));
        memset(    (void*)&partE,0,sizeof(Triplet));
        memset(    (void*)&partP,0,sizeof(Triplet));
        memset((void*)&globPartR,0,sizeof(Triplet));
        memset((void*)&globPartE,0,sizeof(Triplet));
        memset((void*)&globPartP,0,sizeof(Triplet));

        for(int i=0; i!=(int)NeighborCount; ++i)
        {
            NeighborLoad * lo = &nloads[i];
            memset( (void*)lo, 0, sizeof(NeighborLoad));
        }

        sz_nloads=0;
        largest_load_bytes=0;

        break;
    }
    return err;
}

Err_t NEKBONESWMUserCode::conjugateGradient()
{
    Err_t err=0;

    nek_glsc3();

    for(unsigned iter = 0; iter <CGcount; ++iter)
    {
        nek_glsc3();
        nek_gsop("on w");
        nek_glsc3();
        nek_glsc3();
    }

    return err;
}

Err_t NEKBONESWMUserCode::nek_gsop(const char * in_text)
{
    Err_t err=0;

    assert(largest_load_bytes != 0);

    std::vector<char > recvBuffer(largest_load_bytes +1, 0); //+1 for good luck...
    //char * rbuf = &recvBuffer[0];

    SWM_Request emptyRequest = 0;
    std::vector<uint32_t> rrequest(sz_nloads, emptyRequest);

    //In NEKbone, the MPI_Irecv are done first
    for(int i=0; i!=sz_nloads; ++i)
    {
        NeighborLoad b = nloads[i];

        int source = b.rid;

        int tag = 1;  //Each send->recv pair are unique.  So use simplpest tag.

        /*
           mpierr = MPI_Irecv( (void*)rbuf
           , b.load
           , MPI_CHAR
           , source
           , tag
           , MPI_COMM_WORLD
           , rreq
           ); CHECKMPIERR;
           */
        SWM_Irecv(
            source,
            SWM_COMM_WORLD,
            tag,
            /*(BUF_TYPE)rbuf*/ NO_BUFFER,
            &rrequest[i]
        );
    }

    std::vector<char > sendBuffer(largest_load_bytes +1, 0); //+1 for good luck...
    //char * sbuf = &sendBuffer[0];

    uint32_t srequest; //I'm not going to bother with tracking send requests.
    //I hope that is OK.

    //In NEKbone, the MPI_Isend are done second
    for(int i=0; i!=sz_nloads; ++i)
    {
        NeighborLoad b = nloads[i];

        int dest = b.rid;

        int tag = 1; //Each send->recv pair are unique.  So use simplpest tag.

        /*
           mpierr = MPI_Isend( (void*)sbuf
           , b.load
           , MPI_CHAR
           , dest
           , tag
           , MPI_COMM_WORLD
           , &srequest
           ); CHECKMPIERR;
           */
        SWM_Isend(
            dest,
            SWM_COMM_WORLD,
            tag,
            request_vc,
            response_vc,
            /*(BUF_TYPE)sbuf*/ NO_BUFFER,
            b.load,
            pkt_rsp_bytes,
            &srequest,
            req_rt,
            rsp_rt
        );
    }

    {
        //Waiting for all MPI_Recv to be done
        //SWM_Request * rreqall = &rrequest[0];
        int sz = rrequest.size();
        //MPI_Waitall(sz, rreqall, MPI_STATUSES_IGNORE);
        SWM_Waitall(sz, &rrequest[0]);
    }

    return err;
}

Err_t NEKBONESWMUserCode::nek_glsc3()
{
    Err_t err=0;

    //double qx[1];
    //double qres[1];
    //qx[0]=1;
    //qres[0]=0.0;

    //mpierr = MPI_Allreduce( qx, qres, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD );
    //SWM_Allreduce(sizeof(double), 0, comm_world, request_vc, response_vc, (BUF_TYPE)qx, (BUF_TYPE)qres);
    SWM_Allreduce(
        sizeof(double), // payload
        pkt_rsp_bytes, // pkt_rsp_bytes
        SWM_COMM_WORLD,
        request_vc,
        response_vc,
        /*(BUF_TYPE)qx*/ NO_BUFFER,
        /*(BUF_TYPE)qres*/ NO_BUFFER,
        AUTO,
        NULL,
        req_rt,
        rsp_rt
    );

    //double result = qres[0];

    return err;
}

DLL_POSTAMBLE(NEKBONESWMUserCode)
