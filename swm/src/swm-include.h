#ifndef SWM_INCLUDE_H
#define SWM_INCLUDE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define SWM_COMM_WORLD 0
#define NO_BUFFER 0
#define AUTOMATIC 0
#define AUTO 0

typedef int SWM_COMM_ID;
typedef int SWM_PEER;
typedef int SWM_TAG;
typedef int SWM_BUF;
typedef int SWM_ROUTING_TYPE;
typedef uint32_t SWM_BYTES;
typedef int SWM_VC;
typedef int BUF_TYPE;
/*TODO: Find out the right values for SWM_UNKNOWN*/
typedef int SWM_UNKNOWN;
typedef const char* SWM_UNKNOWN2;

struct swm_app_data {
    int final_iteration; // id for final iteration
};

void SWM_Init();

// Call function before first UNION_MPI call that produces any packets. This function will store the information about the application into a space Union can check and respond from.
void SWM_Pass_app_data(struct swm_app_data *);

/*
 * peer: the receiving peer id 
 * comm_id: the communicator id being used
 * tag: tag id 
 * reqvc: virtual channel being used by the message (to be ignored)
 * rspvc: virtual channel being used by the message (to be ignored)
 * buf: the address of sender's buffer in memory
 * bytes: number of bytes to be sent 
 * reqrt and rsprt: routing types (to be ignored) */

void SWM_Send(SWM_PEER peer,
              SWM_COMM_ID comm_id,
              SWM_TAG tag,
              SWM_VC reqvc,
              SWM_VC rspvc,
              SWM_BUF buf,
              SWM_BYTES bytes,
              SWM_BYTES pktrspbytes = 0,
              SWM_ROUTING_TYPE reqrt = AUTOMATIC,
              SWM_ROUTING_TYPE rsprt = AUTOMATIC);

void SWM_Isend(SWM_PEER peer,
              SWM_COMM_ID comm_id,
              SWM_TAG tag,
              SWM_VC reqvc,
              SWM_VC rspvc,
              SWM_BUF buf,
              SWM_BYTES bytes,
              SWM_BYTES pktrspbytes,
              uint32_t * handle,
              SWM_ROUTING_TYPE reqrt = AUTOMATIC,
              SWM_ROUTING_TYPE rsprt = AUTOMATIC);

void SWM_Barrier(
        SWM_COMM_ID comm_id,
        SWM_VC reqvc,
        SWM_VC rspvc,
        SWM_BUF buf, 
        SWM_UNKNOWN auto1,
        SWM_UNKNOWN2 auto2,
        SWM_ROUTING_TYPE reqrt = AUTOMATIC, 
        SWM_ROUTING_TYPE rsprt = AUTOMATIC);

void SWM_Recv(SWM_PEER peer,
        SWM_COMM_ID comm_id,
        SWM_TAG tag,
        SWM_BUF buf);

void SWM_Irecv(SWM_PEER peer,
        SWM_COMM_ID comm_id,
        SWM_TAG tag,
        SWM_BUF buf, 
        uint32_t* handle);

void SWM_Compute(long cycle_count);

void SWM_Wait(uint32_t req_id);

void SWM_Waitall(int len, uint32_t * req_ids);

void SWM_Sendrecv(
         SWM_COMM_ID comm_id,
         SWM_PEER sendpeer,
         SWM_TAG sendtag,
         SWM_VC sendreqvc,
         SWM_VC sendrspvc,
         SWM_BUF sendbuf,
         SWM_BYTES sendbytes,
         SWM_BYTES pktrspbytes,
         SWM_PEER recvpeer,
         SWM_TAG recvtag,
         SWM_BUF recvbuf,
         SWM_ROUTING_TYPE reqrt = AUTOMATIC,
         SWM_ROUTING_TYPE rsprt = AUTOMATIC);

void SWM_Allreduce(
        SWM_BYTES bytes,
        SWM_BYTES rspbytes,
        SWM_COMM_ID comm_id,
        SWM_VC sendreqvc,
        SWM_VC sendrspvc,
        SWM_BUF sendbuf,
        SWM_BUF rcvbuf);

void SWM_Allreduce(
        SWM_BYTES bytes,
        SWM_BYTES rspbytes,
        SWM_COMM_ID comm_id,
        SWM_VC sendreqvc,
        SWM_VC sendrspvc,
        SWM_BUF sendbuf,
        SWM_BUF rcvbuf,
        SWM_UNKNOWN auto1,
        SWM_UNKNOWN2 auto2,
        SWM_ROUTING_TYPE reqrt,
        SWM_ROUTING_TYPE rsprt);

void SWM_Mark_Iteration(
        SWM_TAG iter_tag);

void SWM_Finalize();

#endif
