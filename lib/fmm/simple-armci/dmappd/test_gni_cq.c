/* 
 * CQWrite test 
 *
 * $HeadURL: https://svn.us.cray.com/svn/baker/packages/ugni/trunk/tests/utils/cq_write_test.c $
 * $LastChangedRevision: 1752 $
 */

#include <stdio.h>
#include <stdint.h>
#include <fcntl.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <assert.h>

#ifdef __CRAYXE
#include "gni_pub.h"
#include "pmi.h"

int cookie;
int modes        = 0;

#define TRANS_LEN          128
#define TRANS_LEN_IN_BYTES ((TRANS_LEN)*sizeof(uint64_t))
#define NTRANS             20
#define CACHELINE_MASK     0x3F    /* 64 byte cacheline */

typedef struct {
        gni_mem_handle_t mdh;
        uint64_t addr;
} mdh_addr_t;

unsigned int *MPID_UGNI_AllAddr;

static void
allgather(void *in,void *out, int len)
{
    static int *ivec_ptr=NULL,already_called=0,job_size=0;
    int i,rc;
    int my_rank;
    char *tmp_buf,*out_ptr;

    if(!already_called) {

        rc = PMI_Get_size(&job_size);
        assert(rc == PMI_SUCCESS);
        rc = PMI_Get_rank(&my_rank);
        assert(rc == PMI_SUCCESS);

        ivec_ptr = (int *)malloc(sizeof(int) * job_size);
        assert(ivec_ptr != NULL);

        rc = PMI_Allgather(&my_rank,ivec_ptr,sizeof(int));
        assert(rc == PMI_SUCCESS);

        already_called = 1;

    }

    tmp_buf = (char *)malloc(job_size * len);
    assert(tmp_buf);

    rc = PMI_Allgather(in,tmp_buf,len);
    assert(rc == PMI_SUCCESS);

    out_ptr = out;

    for(i=0;i<job_size;i++) {

        memcpy(&out_ptr[len * ivec_ptr[i]],&tmp_buf[i * len],len);

    }

    free(tmp_buf);
}

static unsigned int
get_gni_nic_address(int device_id)
{
        unsigned int address, cpu_id;
        gni_return_t status;
        int i,alps_dev_id=-1,alps_address=-1;
        char *token,*p_ptr;

        p_ptr = getenv("PMI_GNI_DEV_ID");
        if (!p_ptr) {
                status = GNI_CdmGetNicAddress(device_id, &address, &cpu_id);
                if(status != GNI_RC_SUCCESS) {
                      fprintf(stderr,"FAILED:GNI_CdmGetNicAddress returned error %d\n",status);
                      abort();
                }
        } else {
                while ((token = strtok(p_ptr,":")) != NULL) {
                        alps_dev_id = atoi(token);
                        if (alps_dev_id == device_id) {
                                break;
                        }
                        p_ptr = NULL;
                }
                assert(alps_dev_id != -1);
                p_ptr = getenv("PMI_GNI_LOC_ADDR");
                assert(p_ptr != NULL);
                i = 0;
                while ((token = strtok(p_ptr,":")) != NULL) {
                        if (i == alps_dev_id) {
                                alps_address = atoi(token);
                                break;
                        }
                        p_ptr = NULL;
                        ++i;
                }
                assert(alps_address != -1);
                address = alps_address;
                address = alps_address;
        }

        return address;
}

static void *
gather_nic_addresses(void)
{
    unsigned int local_addr,*alladdrs;
    int size,rc;
    size_t addr_len;

    rc = PMI_Get_size(&size);
    assert(rc == PMI_SUCCESS);

    /*
     * just assume a single gemini device
     */
    local_addr = get_gni_nic_address(0);

    addr_len = sizeof(unsigned int);

    alladdrs = (unsigned int *)malloc(addr_len * size);
    assert(alladdrs != NULL);

    allgather(&local_addr,alladdrs,sizeof(int));

    return (void *)alladdrs;

}

static uint8_t get_ptag(void)
{
        char *p_ptr,*token;
        uint8_t ptag;

        p_ptr = getenv("PMI_GNI_PTAG");
        assert(p_ptr != NULL);  /* something wrong like we haven't called PMI_Init */
        token = strtok(p_ptr,":");
        ptag = (uint8_t)atoi(token);
        return ptag;
}

static uint32_t get_cookie(void)
{
        uint32_t cookie;
        char *p_ptr,*token;

        p_ptr = getenv("PMI_GNI_COOKIE");
        assert(p_ptr != NULL);
        token = strtok(p_ptr,":");
        cookie = (uint32_t)atoi(token);

        return cookie;
}
#endif


int main(int argc,char **argv)
{
#ifdef __CRAYXE
        register int          i;
        int                   rc;
        int                   inst_id, nranks;
        int                   n_entries = 1024;
        int                   device_id = 0;
        int                   send_to, recv_from;
        int                   events_returned;
        int                   test_id;
        unsigned int          local_addr;
        unsigned int          remote_addr;
        gni_cdm_handle_t      cdm_hndl;
        gni_nic_handle_t      nic_hndl;
        gni_return_t          status = GNI_RC_SUCCESS;
        uint32_t              vmdh_index = -1;
        uint64_t              *recv_buffer;
        gni_mem_handle_t      my_remote_mdh;
        uint8_t                  ptag;
        gni_post_descriptor_t fma_data_desc[NTRANS];
        gni_post_descriptor_t *post_desc_ptr;
        mdh_addr_t            my_mdh_addr;
        mdh_addr_t            *remote_mdh_addr_vec;
        gni_cq_handle_t       dst_cq_hndl = NULL;
        gni_cq_handle_t       cq_hndl;
        gni_ep_handle_t       *ep_hndl_array;
        gni_cq_entry_t        event_data;
        int first_spawned;

        if (argc > 1) {
                device_id = atoi(argv[1]);
                fprintf(stderr, "going to try to attach to device %d\n", device_id);
        }

        rc = PMI_Init(&first_spawned);
        assert(rc == PMI_SUCCESS);

        rc = PMI_Get_size(&nranks);
        assert(rc == PMI_SUCCESS);

        rc = PMI_Get_rank(&inst_id);
        assert(rc == PMI_SUCCESS);

        ptag = get_ptag();
        cookie = get_cookie();

        /* Create and attach to the communication domain. */

        status = GNI_CdmCreate(inst_id, ptag, cookie, modes, &cdm_hndl);
        if (status != GNI_RC_SUCCESS) {
                fprintf(stderr, "FAIL: GNI_CdmCreate returned error %s\n", gni_err_str[status]);
                PMI_Abort(-1,"pmi abort called");
        } else {
                fprintf(stderr, "Inst_id %d created CDM, now attaching...\n", inst_id);
        }

        status = GNI_CdmAttach(cdm_hndl, device_id, &local_addr, &nic_hndl);
        if (status != GNI_RC_SUCCESS) {
                fprintf(stderr, "FAIL: GNI_CdmAttach returned error %s\n", gni_err_str[status]);
                PMI_Abort(-1,"pmi abort called");
        } else {
                fprintf(stderr, "Inst_id %d attached CDM to NIC\n", inst_id);
        }

        /* Create the local completion queue */

        status = GNI_CqCreate(nic_hndl, n_entries, 0, GNI_CQ_NOBLOCK, NULL, NULL, &cq_hndl);
        if (status != GNI_RC_SUCCESS) {
                fprintf(stderr, "FAIL: GNI_CqCreate returned error %s\n", gni_err_str[status]);
                PMI_Abort(-1,"pmi abort called");
        } else {
                fprintf(stderr, "Inst_id %d created CQ\n", inst_id);
        }


        /* Create the destination completion queue for receiving micro-messages,
           make this queue considerably larger than number of transfers */

        status = GNI_CqCreate(nic_hndl, NTRANS*10, 0, GNI_CQ_NOBLOCK, NULL, NULL, &dst_cq_hndl);
        if (status != GNI_RC_SUCCESS) {
                fprintf(stderr, "FAIL: GNI_CqCreate returned error %s\n", gni_err_str[status]);
                PMI_Abort(-1,"pmi abort called");
        } else {
                fprintf(stderr, "Inst_id %d created CQ\n", inst_id);
        }

        /* Create the endpoints. They need to be bound to allow later CQwrites to 
           them. */

        ep_hndl_array = (gni_ep_handle_t *)malloc(nranks * sizeof(gni_ep_handle_t));
        assert(ep_hndl_array != NULL);

        MPID_UGNI_AllAddr = (unsigned int *)gather_nic_addresses();

        for (i=0; i<nranks; i++) {
                status = GNI_EpCreate(nic_hndl, cq_hndl, &ep_hndl_array[i]);
                if(status != GNI_RC_SUCCESS) {
                        fprintf(stderr, "FAIL: GNI_EpCreate returned error %s\n", gni_err_str[status]);
                        PMI_Abort(-1,"pmi abort called");
                }
                remote_addr = MPID_UGNI_AllAddr[i];
                status = GNI_EpBind(ep_hndl_array[i], remote_addr, i);
                if(status != GNI_RC_SUCCESS) {
                        fprintf(stderr, "FAIL: GNI_EpCreate returned error %s\n", gni_err_str[status]);
                        PMI_Abort(-1,"pmi abort called");
                }
        }

        /* Allocate a 4096 receive buffer to hook the destination cq to */

        recv_buffer = (uint64_t *)calloc(512, sizeof(uint64_t));
        assert(recv_buffer != NULL);

        status = GNI_MemRegister(nic_hndl, (uint64_t)recv_buffer,
                                 4096, dst_cq_hndl, 
                                 GNI_MEM_READWRITE | GNI_MEM_USE_GART, vmdh_index,
                                 &my_remote_mdh);
        if (status != GNI_RC_SUCCESS) {
                fprintf(stderr,"GNI_MemRegister returned error %s\n", gni_err_str[status]);
                PMI_Abort(-1,"pmi abort called");
        }

        /* Gather up all of the mdh's over the socket network, 
           this also serves as a barrier */

        remote_mdh_addr_vec = (mdh_addr_t *)malloc(nranks * sizeof(mdh_addr_t));
        assert(remote_mdh_addr_vec);

        my_mdh_addr.addr = (uint64_t)recv_buffer;
        my_mdh_addr.mdh = my_remote_mdh;

        allgather(&my_mdh_addr, remote_mdh_addr_vec,sizeof(mdh_addr_t));

        send_to = (inst_id + 1)% nranks;
        recv_from = (nranks + inst_id - 1)% nranks;

        fprintf(stderr, "Rank %d sending destination CQ writes\n", inst_id);
        for (i=0; i<NTRANS; i++) {

                /* Prepare and send the data, we send with hashing adaption
                   disabled so that cq data arrives in order at target node*/

                fma_data_desc[i].type = GNI_POST_CQWRITE;
                fma_data_desc[i].cq_mode = GNI_CQMODE_GLOBAL_EVENT;
                fma_data_desc[i].dlvr_mode = GNI_DLVMODE_NO_ADAPT;
                fma_data_desc[i].cqwrite_value = 0xdeadbeefdeadbeef + i;
                fma_data_desc[i].remote_mem_hndl = remote_mdh_addr_vec[send_to].mdh;

                status = GNI_PostCqWrite(ep_hndl_array[send_to], &fma_data_desc[i]);
                if (status != GNI_RC_SUCCESS) {
                        fprintf(stderr, "FAIL: GNI_PostCqWrite returned error %d\n", status);
                        PMI_Abort(-1,"pmi abort called");
                } 

                status = GNI_RC_NOT_DONE;
                while (status == GNI_RC_NOT_DONE) {

                        status = GNI_CqGetEvent(dst_cq_hndl, &event_data);
                        if (status == GNI_RC_SUCCESS) {
                            if (GNI_CQ_OVERRUN(event_data)) {
                                    fprintf(stderr, "FAIL: erroneous CQ_OVERRUN detected\n");
                                    PMI_Abort(-1,"pmi abort called");
                            }
                            if (GNI_CQ_GET_DATA(event_data) !=
                                   GNI_CQ_GET_DATA(fma_data_desc[i].cqwrite_value)) {
                                        fprintf(stderr, "FAIL: erroneous CQ value detected\n");
                                        PMI_Abort(-1,"pmi abort called");
                            }
                        }
                }
                        
        } /* end loop over NTRANS */

        /*
         * now get all of my events back
         */
        events_returned = 0;
        fprintf(stderr, "Rank %d data transfers complete, checking local CQ events\n",
                inst_id);

        while (events_returned < NTRANS) {
                status = GNI_CqGetEvent(cq_hndl, &event_data);
                if (status == GNI_RC_SUCCESS) {
                        if (GNI_CQ_OVERRUN(event_data)) {
                                fprintf(stderr, "FAIL: erroneous CQ_OVERRUN detected\n");
                                PMI_Abort(-1,"pmi abort called");
                        }

                        test_id = GNI_CQ_GET_INST_ID(event_data);
                        if (test_id != send_to) {
                                fprintf(stderr, "FAIL: inst_id %d got incorrect inst_id(%d) in event_data\n",
                                        inst_id,test_id);
                                PMI_Abort(-1,"pmi abort called");
                        }

                        status = GNI_GetCompleted(cq_hndl, event_data, &post_desc_ptr);
                        if (status != GNI_RC_SUCCESS) {
                                fprintf(stderr, "FAIL: got error return from GNI_Get_completed %s\n", 
                                        gni_err_str[status]);
                                PMI_Abort(-1,"pmi abort called");
                        }
                        events_returned++;
                }
                else if (status != GNI_RC_NOT_DONE) {
                        fprintf(stderr, "FAIL: error return from GNI_CqGetEvent %s\n",
                                gni_err_str[status]);
                        PMI_Abort(-1,"pmi abort called");
                }
                else {  /* data has not been sent yet; do nothing */
                }
        }

        fprintf(stderr, "Rank %d CQ events received, test passed\n", inst_id);

        PMI_Finalize();
#endif
        return 0;
}


