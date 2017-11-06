/* ---------------------------------------------------------------- */
/* (C)Copyright IBM Corp.  2007, 2008                               */
/* ---------------------------------------------------------------- */
/**
 * \file ad_bg_aggrs.c
 * \brief The externally used function from this file is is declared in ad_bg_aggrs.h
 */

/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *   Copyright (C) 1997-2001 University of Chicago.
 *   See COPYRIGHT notice in top-level directory.
 */

/*#define TRACE_ON */

#include "adio.h"
#include "adio_cb_config_list.h"
#include "ad_bg.h"
#include "ad_bg_pset.h"
#include "ad_bg_aggrs.h"
#ifdef AGGREGATION_PROFILE
#include "mpe.h"
#endif

#include "mpidi_macros.h"

#ifdef USE_DBG_LOGGING
  #define AGG_DEBUG 1
#endif

static int aggrsInPsetSize=0;
static int *aggrsInPset=NULL;

/* Comments copied from common:
 * This file contains four functions:
 *
 * ADIOI_Calc_aggregator()
 * ADIOI_Calc_file_domains()
 * ADIOI_Calc_my_req()
 * ADIOI_Calc_others_req()
 *
 * The last three of these were originally in ad_read_coll.c, but they are
 * also shared with ad_write_coll.c.  I felt that they were better kept with
 * the rest of the shared aggregation code.
 */

/* Discussion of values available from above:
 *
 * ADIO_Offset st_offsets[0..nprocs-1]
 * ADIO_Offset end_offsets[0..nprocs-1]
 *    These contain a list of start and end offsets for each process in
 *    the communicator.  For example, an access at loc 10, size 10 would
 *    have a start offset of 10 and end offset of 19.
 * int nprocs
 *    number of processors in the collective I/O communicator
 * ADIO_Offset min_st_offset
 * ADIO_Offset fd_start[0..nprocs_for_coll-1]
 *    starting location of "file domain"; region that a given process will
 *    perform aggregation for (i.e. actually do I/O)
 * ADIO_Offset fd_end[0..nprocs_for_coll-1]
 *    start + size - 1 roughly, but it can be less, or 0, in the case of
 *    uneven distributions
 */

/* forward declaration */
static void
ADIOI_BG_compute_agg_ranklist_serial ( ADIO_File fd,
					const ADIOI_BG_ConfInfo_t *confInfo,
					ADIOI_BG_ProcInfo_t *all_procInfo,
					int *aggrsInPset );

/*
 * Compute the aggregator-related parameters that are required in 2-phase collective IO of ADIO.
 * The parameters are
 * 	. the number of aggregators (proxies) : fd->hints->cb_nodes
 *	. the ranks of the aggregators :        fd->hints->ranklist
 * By compute these two parameters in a BG-PSET-aware way, the default 2-phase collective IO of
 *	ADIO can work more efficiently.
 */
int
ADIOI_BG_gen_agg_ranklist(ADIO_File fd, int n_aggrs_per_pset)
{
    int r, s;
    ADIOI_BG_ProcInfo_t  *procInfo, *all_procInfo;
    ADIOI_BG_ConfInfo_t  *confInfo;
    TRACE_ERR("Entering ADIOI_BG_gen_agg_ranklist\n");

    MPI_Comm_size( fd->comm, &s );
    MPI_Comm_rank( fd->comm, &r );

  /* Collect individual BG personality information */
    confInfo = ADIOI_BG_ConfInfo_new ();
    procInfo = ADIOI_BG_ProcInfo_new ();
    ADIOI_BG_persInfo_init( confInfo, procInfo, s, r, n_aggrs_per_pset, fd->comm);

  /* Gather BG personality infomation onto process 0 */
    /* if (r == 0) */
    all_procInfo  = ADIOI_BG_ProcInfo_new_n  (s);
    if(s > aggrsInPsetSize)
    {
      if(aggrsInPset) ADIOI_Free(aggrsInPset);
      aggrsInPset   = (int *) ADIOI_Malloc (s *sizeof(int));
      aggrsInPsetSize = s;
    }


    MPI_Gather( (void *)procInfo,     sizeof(ADIOI_BG_ProcInfo_t), MPI_BYTE,
		(void *)all_procInfo, sizeof(ADIOI_BG_ProcInfo_t), MPI_BYTE,
		0,
		fd->comm );

  /* Compute a list of the ranks of chosen IO proxy CN on process 0 */
    if (r == 0) {
	ADIOI_BG_compute_agg_ranklist_serial (fd, confInfo, all_procInfo, aggrsInPset);
	/* ADIOI_BG_ProcInfo_free (all_procInfo);*/
    }
    ADIOI_BG_ProcInfo_free (all_procInfo);

  /* Send the info of IO proxy CN to all processes and keep the info in fd->hints struct.
     Declared in adio_cb_config_list.h */
    ADIOI_cb_bcast_rank_map(fd);

  /* Broadcast the BG-GPFS related file domain info */
    MPI_Bcast( (void *)aggrsInPset,
		fd->hints->cb_nodes * sizeof(int), MPI_BYTE,
		0,
		fd->comm );

    ADIOI_BG_persInfo_free( confInfo, procInfo );
    TRACE_ERR("Leaving ADIOI_BG_gen_agg_ranklist\n");
    return 0;
}


/* There are some number of bridge nodes (randomly) distributed through the job
 * We need to split the nodes among the bridge nodes */
/* Maybe find which bridge node is closer (manhattan distance) and try to
 * distribute evenly.
 */
/*
 * Pick IO aggregators based on the under PSET organization and stores the ranks of the proxy CNs in tmp_ranklist.
 * The first order of tmp_ranklist is : PSET number
 * The secondary order of the list is determined in ADIOI_BG_select_agg_in_pset() and thus adjustable.
 */
typedef struct
{
   int rank;
   int bridge;
} sortstruct;

static int intsort(const void *p1, const void *p2)
{
   sortstruct *i1, *i2;
   i1 = (sortstruct *)p1;
   i2 = (sortstruct *)p2;
   return(i1->bridge - i2->bridge);
}

static int
ADIOI_BG_compute_agg_ranklist_serial_do (const ADIOI_BG_ConfInfo_t *confInfo,
					  ADIOI_BG_ProcInfo_t       *all_procInfo,
					  int *aggrsInPset,
					  int *tmp_ranklist)
{
    TRACE_ERR("Entering ADIOI_BG_compute_agg_ranklist_serial_do\n");
   /* BES: This should be done in the init routines probably. */
    int i, j;
    int aggTotal;
    int distance, numAggs;
    int *aggList;
    /* Aggregators will be midpoints between sorted MPI rank lists of who shares a given
     * bridge node */

   sortstruct *bridgelist = (sortstruct *)ADIOI_Malloc(confInfo->nProcs * sizeof(sortstruct));
   for(i=0; i < confInfo->nProcs; i++)
   {
      bridgelist[i].bridge = all_procInfo[i].bridgeRank;
      bridgelist[i].rank = i;
      TRACE_ERR("bridgelist[%d].bridge: %d .rank: %d\n", i, bridgelist[i].bridge, i);
   }

   /* This list contains rank->bridge info. Now, we need to sort this list. */
   qsort(bridgelist, confInfo->nProcs, sizeof(sortstruct), intsort);

   /* In this array, we can pick an appropriate number of midpoints based on
    * our bridgenode index and the number of aggregators */

   numAggs = confInfo->aggRatio * confInfo->ioMaxSize /*virtualPsetSize*/;
   if(numAggs == 1)
      aggTotal = 1;
   else
   /* the number of aggregators is (numAggs per bridgenode) plus each
    * bridge node is an aggregator */
      aggTotal = confInfo->numBridgeRanks * (numAggs+1);

   distance = (confInfo->ioMaxSize /*virtualPsetSize*/ / numAggs);
   TRACE_ERR("numBridgeRanks: %d, aggRatio: %f numBridge: %d pset size: %d numAggs: %d distance: %d, aggTotal: %d\n", confInfo->numBridgeRanks, confInfo->aggRatio, confInfo->numBridgeRanks,  confInfo->ioMaxSize /*virtualPsetSize*/, numAggs, distance, aggTotal);
   aggList = (int *)ADIOI_Malloc(aggTotal * sizeof(int));


   /* For each bridge node, determine who the aggregators will be */
   /* basically, the n*distance and bridge node */
   if(aggTotal == 1) /* special case when we only have one bridge node */
      aggList[0] = bridgelist[0].bridge;
   else
   {
      for(i=0; i < confInfo->numBridgeRanks; i++)
      {
         aggList[i]=bridgelist[i*confInfo->ioMaxSize /*virtualPsetSize*/].bridge;
         TRACE_ERR("aggList[%d]: %d\n", i, aggList[i]);

         for(j = 0; j < numAggs; j++)
         {
            /* Sets up a list of nodes which will act as aggregators. numAggs
             * per bridge node total. The list of aggregators is
             * bridgeNodes
             * bridgeNode[0]aggr[0]
             * bridgeNode[0]aggr[1]...
             * bridgeNode[0]aggr[N]...
             * ...
             * bridgeNode[N]aggr[0]..
             * bridgeNode[N]aggr[N]
             */
            aggList[i*numAggs+j+confInfo->numBridgeRanks] = bridgelist[i*confInfo->ioMaxSize /*virtualPsetSize*/ + j*distance+1].rank;
            TRACE_ERR("(post bridge) agglist[%d] -> %d\n", confInfo->numBridgeRanks +i*numAggs+j, aggList[i*numAggs+j+confInfo->numBridgeRanks]);
         }
      }
   }

   memcpy(tmp_ranklist, aggList, (numAggs*confInfo->numBridgeRanks+numAggs)*sizeof(int));
   for(i=0;i<aggTotal;i++)
   {
      TRACE_ERR("tmp_ranklist[%d]: %d\n", i, tmp_ranklist[i]);
   }


   ADIOI_Free (bridgelist);
   ADIOI_Free (aggList);

   TRACE_ERR("Leaving ADIOI_BG_compute_agg_ranklist_serial_do\n");
   return aggTotal;

}

/*
 * compute aggregators ranklist and put it into fd->hints struct
 */
static void
ADIOI_BG_compute_agg_ranklist_serial ( ADIO_File fd,
					const ADIOI_BG_ConfInfo_t *confInfo,
					ADIOI_BG_ProcInfo_t *all_procInfo,
					int *aggrsInPset )
{
    TRACE_ERR("Entering ADIOI_BG_compute_agg_ranklist_serial\n");
    int i;
    int naggs;
    int size;
    int *tmp_ranklist;

  /* compute the ranklist of IO aggregators and put into tmp_ranklist */
    tmp_ranklist = (int *) ADIOI_Malloc (confInfo->nProcs * sizeof(int));

#   if AGG_DEBUG
    for (i=0; i<confInfo->nProcs; i++) {
      DBG_FPRINTF(stderr, "\tcpuid %1d, rank = %6d\n", all_procInfo[i].coreID, all_procInfo[i].rank );
    }
#   endif

    naggs=
    ADIOI_BG_compute_agg_ranklist_serial_do (confInfo, all_procInfo, aggrsInPset, tmp_ranklist);

#   define VERIFY 1
#   if VERIFY
    DBG_FPRINTF(stderr, "\tconfInfo = min: %3d, max: %3d, naggrs: %3d, bridge: %3d, nprocs: %3d, vpset: %3d, tsize: %3d, ratio: %.4f; naggs = %d\n",
	    confInfo->ioMinSize        ,
	    confInfo->ioMaxSize        ,
	    confInfo->nAggrs           ,
	    confInfo->numBridgeRanks ,
	    confInfo->nProcs          ,
	    confInfo->ioMaxSize /*virtualPsetSize*/          ,
      confInfo->cpuIDsize,
	    confInfo->aggRatio        ,
	    naggs );
#   endif
    MPI_Comm_size( fd->comm, &size );
    /* This fix is for when the bridgenode rnk is not part of the particular
     * subcomm associated with this MPI File operation. I don't know if
     * this is the best/right answer but it passes the test cases at least.
     * I don't know how common file IO in subcomms is anyway... */
    for(i=0;i<naggs;i++)
    {
      if(tmp_ranklist[i] > size)
      {
         TRACE_ERR("Using 0 as tmp_ranklist[%d] instead of %d for comm %x\n",
               i, tmp_ranklist[i], fd->comm);
         tmp_ranklist[i] = 0;
      }
   }

#   if AGG_DEBUG
    for (i=0; i<naggs; i++) {
      DBG_FPRINTF(stderr, "\taggr %-4d = %6d\n", i, tmp_ranklist[i] );
    }
#   endif

  /* copy the ranklist of IO aggregators to fd->hints */
    if(fd->hints->ranklist != NULL) ADIOI_Free (fd->hints->ranklist);

    fd->hints->cb_nodes = naggs;
    fd->hints->ranklist = (int *) ADIOI_Malloc (naggs * sizeof(int));
    memcpy( fd->hints->ranklist, tmp_ranklist, naggs*sizeof(int) );

  /* */
    ADIOI_Free( tmp_ranklist );
    TRACE_ERR("Leaving ADIOI_BG_compute_agg_ranklist_serial\n");
    return;
}

/* Description from common/ad_aggregate.c.  (Does it completely apply to bg?)
 * ADIOI_Calc_aggregator()
 *
 * The intention here is to implement a function which provides basically
 * the same functionality as in Rajeev's original version of
 * ADIOI_Calc_my_req().  He used a ceiling division approach to assign the
 * file domains, and we use the same approach here when calculating the
 * location of an offset/len in a specific file domain.  Further we assume
 * this same distribution when calculating the rank_index, which is later
 *  used to map to a specific process rank in charge of the file domain.
 *
 * A better (i.e. more general) approach would be to use the list of file
 * domains only.  This would be slower in the case where the
 * original ceiling division was used, but it would allow for arbitrary
 * distributions of regions to aggregators.  We'd need to know the
 * nprocs_for_coll in that case though, which we don't have now.
 *
 * Note a significant difference between this function and Rajeev's old code:
 * this code doesn't necessarily return a rank in the range
 * 0..nprocs_for_coll; instead you get something in 0..nprocs.  This is a
 * result of the rank mapping; any set of ranks in the communicator could be
 * used now.
 *
 * Returns an integer representing a rank in the collective I/O communicator.
 *
 * The "len" parameter is also modified to indicate the amount of data
 * actually available in this file domain.
 */
/*
 * This is more general aggregator search function which does not base on the assumption
 * that each aggregator hosts the file domain with the same size
 */
int ADIOI_BG_Calc_aggregator(ADIO_File fd,
			      ADIO_Offset off,
			      ADIO_Offset min_off,
			      ADIO_Offset *len,
			      ADIO_Offset fd_size,
			      ADIO_Offset *fd_start,
			      ADIO_Offset *fd_end)
{
    int rank_index, rank;
    ADIO_Offset avail_bytes;
    TRACE_ERR("Entering ADIOI_BG_Calc_aggregator\n");

    ADIOI_BG_assert ( (off <= fd_end[fd->hints->cb_nodes-1] && off >= min_off && fd_start[0] >= min_off ) );

    /* binary search --> rank_index is returned */
    int ub = fd->hints->cb_nodes;
    int lb = 0;
    /* get an index into our array of aggregators */
    /* Common code for striping - bg doesn't use it but it's
       here to make diff'ing easier.
    rank_index = (int) ((off - min_off + fd_size)/ fd_size - 1);

    if (fd->hints->striping_unit > 0) {
        * wkliao: implementation for file domain alignment
           fd_start[] and fd_end[] have been aligned with file lock
	   boundaries when returned from ADIOI_Calc_file_domains() so cannot
	   just use simple arithmatic as above *
        rank_index = 0;
        while (off > fd_end[rank_index]) rank_index++;
    }
    bg does it's own striping below
    */
    rank_index = fd->hints->cb_nodes / 2;
    while ( off < fd_start[rank_index] || off > fd_end[rank_index] ) {
	if ( off > fd_end  [rank_index] ) {
	    lb = rank_index;
	    rank_index = (rank_index + ub) / 2;
	}
	else
	if ( off < fd_start[rank_index] ) {
	    ub = rank_index;
	    rank_index = (rank_index + lb) / 2;
	}
    }
    /* we index into fd_end with rank_index, and fd_end was allocated to be no
     * bigger than fd->hins->cb_nodes.   If we ever violate that, we're
     * overrunning arrays.  Obviously, we should never ever hit this abort */
    if (rank_index >= fd->hints->cb_nodes || rank_index < 0) {
        FPRINTF(stderr, "Error in ADIOI_Calc_aggregator(): rank_index(%d) >= fd->hints->cb_nodes (%d) fd_size=%lld off=%lld\n",
			rank_index,fd->hints->cb_nodes,fd_size,off);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    /* DBG_FPRINTF ("ADIOI_BG_Calc_aggregator: rank_index = %d\n",
       rank_index ); */

    /*
     * remember here that even in Rajeev's original code it was the case that
     * different aggregators could end up with different amounts of data to
     * aggregate.  here we use fd_end[] to make sure that we know how much
     * data this aggregator is working with.
     *
     * the +1 is to take into account the end vs. length issue.
     */
    avail_bytes = fd_end[rank_index] + 1 - off;
    if (avail_bytes < *len && avail_bytes > 0) {
        /* this file domain only has part of the requested contig. region */

        *len = avail_bytes;
    }

    /* map our index to a rank */
    /* NOTE: FOR NOW WE DON'T HAVE A MAPPING...JUST DO 0..NPROCS_FOR_COLL */
    rank = fd->hints->ranklist[rank_index];
    TRACE_ERR("Leaving ADIOI_BG_Calc_aggregator\n");

    return rank;
}

/*
 * Compute a dynamic access range based file domain partition among I/O aggregators,
 * which align to the GPFS block size
 * Divide the I/O workload among "nprocs_for_coll" processes. This is
 * done by (logically) dividing the file into file domains (FDs); each
 * process may directly access only its own file domain.
 * Additional effort is to make sure that each I/O aggregator get
 * a file domain that aligns to the GPFS block size.  So, there will
 * not be any false sharing of GPFS file blocks among multiple I/O nodes.
 *
 * The common version of this now accepts a min_fd_size and striping_unit.
 * It doesn't seem necessary here (using GPFS block sizes) but keep it in mind
 * (e.g. we could pass striping unit instead of using fs_ptr->blksize).
 */
void ADIOI_BG_GPFS_Calc_file_domains(ADIO_Offset *st_offsets,
                                      ADIO_Offset *end_offsets,
                                      int          nprocs,
                                      int          nprocs_for_coll,
                                      ADIO_Offset *min_st_offset_ptr,
                                      ADIO_Offset **fd_start_ptr,
                                      ADIO_Offset **fd_end_ptr,
                                      ADIO_Offset *fd_size_ptr,
                                      void        *fs_ptr)
{
    ADIO_Offset min_st_offset, max_end_offset, *fd_start, *fd_end, *fd_size;
    int i, aggr;
    TRACE_ERR("Entering ADIOI_BG_GPFS_Calc_file_domains\n");

#ifdef AGGREGATION_PROFILE
    MPE_Log_event (5004, 0, NULL);
#endif

#   if AGG_DEBUG
    static char myname[] = "ADIOI_BG_GPFS_Calc_file_domains";
    DBG_FPRINTF(stderr, "%s(%d): %d aggregator(s)\n",
	    myname,__LINE__,nprocs_for_coll);
#   endif
    __blksize_t blksize = 1048576; /* default to 1M */
    if(fs_ptr && ((ADIOI_BG_fs*)fs_ptr)->blksize) /* ignore null ptr or 0 blksize */
      blksize = ((ADIOI_BG_fs*)fs_ptr)->blksize;
#   if AGG_DEBUG
    DBG_FPRINTF(stderr,"%s(%d): Blocksize=%ld\n",myname,__LINE__,blksize);
#   endif
/* find min of start offsets and max of end offsets of all processes */
    min_st_offset  = st_offsets [0];
    max_end_offset = end_offsets[0];
    for (i=1; i<nprocs; i++) {
        min_st_offset = ADIOI_MIN(min_st_offset, st_offsets[i]);
        max_end_offset = ADIOI_MAX(max_end_offset, end_offsets[i]);
    }

    /* DBG_FPRINTF(stderr, "_calc_file_domains, min_st_offset, max_
       = %qd, %qd\n", min_st_offset, max_end_offset );*/

    /* determine the "file domain (FD)" of each process, i.e., the portion of
       the file that will be "owned" by each process */

    ADIO_Offset gpfs_ub       = (max_end_offset +blksize-1) / blksize * blksize - 1;
    ADIO_Offset gpfs_lb       = min_st_offset / blksize * blksize;
    ADIO_Offset gpfs_ub_rdoff = (max_end_offset +blksize-1) / blksize * blksize - 1 - max_end_offset;
    ADIO_Offset gpfs_lb_rdoff = min_st_offset - min_st_offset / blksize * blksize;
    ADIO_Offset fd_gpfs_range = gpfs_ub - gpfs_lb + 1;

    int         naggs    = nprocs_for_coll;

    /* Tweak the file domains so that no fd is smaller than a threshold.  We
     * have to strike a balance between efficency and parallelism: somewhere
     * between 10k processes sending 32-byte requests and one process sending a
     * 320k request is a (system-dependent) sweet spot

    This is from the common code - the new min_fd_size parm that we didn't implement.
    (And common code uses a different declaration of fd_size so beware)

    if (fd_size < min_fd_size)
        fd_size = min_fd_size;
    */
    fd_size              = (ADIO_Offset *) ADIOI_Malloc(nprocs_for_coll * sizeof(ADIO_Offset));
    *fd_start_ptr        = (ADIO_Offset *) ADIOI_Malloc(nprocs_for_coll * sizeof(ADIO_Offset));
    *fd_end_ptr          = (ADIO_Offset *) ADIOI_Malloc(nprocs_for_coll * sizeof(ADIO_Offset));
    fd_start             = *fd_start_ptr;
    fd_end               = *fd_end_ptr;

    ADIO_Offset n_gpfs_blk    = fd_gpfs_range / blksize;
    ADIO_Offset nb_cn_small   = n_gpfs_blk/naggs;
    ADIO_Offset naggs_large   = n_gpfs_blk - naggs * (n_gpfs_blk/naggs);
    ADIO_Offset naggs_small   = naggs - naggs_large;

    for (i=0; i<naggs; i++)
        if (i < naggs_small) fd_size[i] = nb_cn_small     * blksize;
                        else fd_size[i] = (nb_cn_small+1) * blksize;

#   if AGG_DEBUG
     DBG_FPRINTF(stderr,"%s(%d): "
                   "gpfs_ub       %llu, "
                   "gpfs_lb       %llu, "
                   "gpfs_ub_rdoff %llu, "
                   "gpfs_lb_rdoff %llu, "
                   "fd_gpfs_range %llu, "
                   "n_gpfs_blk    %llu, "
                   "nb_cn_small   %llu, "
                   "naggs_large   %llu, "
                   "naggs_small   %llu, "
                   "\n",
                   myname,__LINE__,
                   gpfs_ub      ,
                   gpfs_lb      ,
                   gpfs_ub_rdoff,
                   gpfs_lb_rdoff,
                   fd_gpfs_range,
                   n_gpfs_blk   ,
                   nb_cn_small  ,
                   naggs_large  ,
                   naggs_small
                   );
#   endif

    fd_size[0]       -= gpfs_lb_rdoff;
    fd_size[naggs-1] -= gpfs_ub_rdoff;

    /* compute the file domain for each aggr */
    ADIO_Offset offset = min_st_offset;
    for (aggr=0; aggr<naggs; aggr++) {
        fd_start[aggr] = offset;
        fd_end  [aggr] = offset + fd_size[aggr] - 1;
        offset += fd_size[aggr];
    }

    *fd_size_ptr = fd_size[0];
    *min_st_offset_ptr = min_st_offset;

#ifdef AGGREGATION_PROFILE
    MPE_Log_event (5005, 0, NULL);
#endif
    ADIOI_Free (fd_size);
    TRACE_ERR("Leaving ADIOI_BG_GPFS_Calc_file_domains\n");
}

/*
 * When a process is an IO aggregator, this will return its index in the aggrs list.
 * Otherwise, this will return -1
 */
int ADIOI_BG_Aggrs_index( ADIO_File fd, int myrank )
{
    int i;
    for (i=0; i<fd->hints->cb_nodes; i++)
	if (fd->hints->ranklist[i] == myrank) return i;
    return -1;
}

/*
 * ADIOI_BG_Calc_my_req() overrides ADIOI_Calc_my_req for the default implementation
 * is specific for static file domain partitioning.
 *
 * ADIOI_Calc_my_req() - calculate what portions of the access requests
 * of this process are located in the file domains of various processes
 * (including this one)
 */
void ADIOI_BG_Calc_my_req(ADIO_File fd, ADIO_Offset *offset_list, ADIO_Offset *len_list,
			   int contig_access_count, ADIO_Offset
			   min_st_offset, ADIO_Offset *fd_start,
			   ADIO_Offset *fd_end, ADIO_Offset fd_size,
			   int nprocs,
			   int *count_my_req_procs_ptr,
			   int **count_my_req_per_proc_ptr,
			   ADIOI_Access **my_req_ptr,
			   int **buf_idx_ptr)
/* Possibly reconsider if buf_idx's are ok as int's, or should they be aints/offsets?
   They are used as memory buffer indices so it seems like the 2G limit is in effect */
{
    int *count_my_req_per_proc, count_my_req_procs, *buf_idx;
    int i, l, proc;
    ADIO_Offset fd_len, rem_len, curr_idx, off;
    ADIOI_Access *my_req;
    TRACE_ERR("Entering ADIOI_BG_Calc_my_req\n");

#ifdef AGGREGATION_PROFILE
    MPE_Log_event (5024, 0, NULL);
#endif

    *count_my_req_per_proc_ptr = (int *) ADIOI_Calloc(nprocs,sizeof(int));
    count_my_req_per_proc = *count_my_req_per_proc_ptr;
/* count_my_req_per_proc[i] gives the no. of contig. requests of this
   process in process i's file domain. calloc initializes to zero.
   I'm allocating memory of size nprocs, so that I can do an
   MPI_Alltoall later on.*/

    buf_idx = (int *) ADIOI_Malloc(nprocs*sizeof(int));
/* buf_idx is relevant only if buftype_is_contig.
   buf_idx[i] gives the index into user_buf where data received
   from proc. i should be placed. This allows receives to be done
   without extra buffer. This can't be done if buftype is not contig. */

    /* initialize buf_idx to -1 */
    for (i=0; i < nprocs; i++) buf_idx[i] = -1;

    /* one pass just to calculate how much space to allocate for my_req;
     * contig_access_count was calculated way back in ADIOI_Calc_my_off_len()
     */
    for (i=0; i < contig_access_count; i++) {
	/* short circuit offset/len processing if len == 0
	 * 	(zero-byte  read/write */
	if (len_list[i] == 0)
		continue;
	off = offset_list[i];
	fd_len = len_list[i];
	/* note: we set fd_len to be the total size of the access.  then
	 * ADIOI_Calc_aggregator() will modify the value to return the
	 * amount that was available from the file domain that holds the
	 * first part of the access.
	 */
  /* BES */
	proc = ADIOI_BG_Calc_aggregator(fd, off, min_st_offset, &fd_len, fd_size,
				     fd_start, fd_end);
	count_my_req_per_proc[proc]++;

	/* figure out how much data is remaining in the access (i.e. wasn't
	 * part of the file domain that had the starting byte); we'll take
	 * care of this data (if there is any) in the while loop below.
	 */
	rem_len = len_list[i] - fd_len;

	while (rem_len > 0) {
	    off += fd_len; /* point to first remaining byte */
	    fd_len = rem_len; /* save remaining size, pass to calc */
	    proc = ADIOI_BG_Calc_aggregator(fd, off, min_st_offset, &fd_len,
					 fd_size, fd_start, fd_end);

	    count_my_req_per_proc[proc]++;
	    rem_len -= fd_len; /* reduce remaining length by amount from fd */
	}
    }

/* now allocate space for my_req, offset, and len */

    *my_req_ptr = (ADIOI_Access *)
	ADIOI_Malloc(nprocs*sizeof(ADIOI_Access));
    my_req = *my_req_ptr;

    count_my_req_procs = 0;
    for (i=0; i < nprocs; i++) {
	if (count_my_req_per_proc[i]) {
	    my_req[i].offsets = (ADIO_Offset *)
		ADIOI_Malloc(count_my_req_per_proc[i] * sizeof(ADIO_Offset));
	    my_req[i].lens = (int *)
		ADIOI_Malloc(count_my_req_per_proc[i] * sizeof(int));
	    count_my_req_procs++;
	}
	my_req[i].count = 0;  /* will be incremented where needed
				      later */
    }

/* now fill in my_req */
    curr_idx = 0;
    for (i=0; i<contig_access_count; i++) {
	/* short circuit offset/len processing if len == 0
	 * 	(zero-byte  read/write */
	if (len_list[i] == 0)
		continue;
	off = offset_list[i];
	fd_len = len_list[i];
	proc = ADIOI_BG_Calc_aggregator(fd, off, min_st_offset, &fd_len, fd_size,
				     fd_start, fd_end);

	/* for each separate contiguous access from this process */
	if (buf_idx[proc] == -1)
  {
    ADIOI_Assert(curr_idx == (int) curr_idx);
    buf_idx[proc] = (int) curr_idx;
  }

	l = my_req[proc].count;
	curr_idx += fd_len;

	rem_len = len_list[i] - fd_len;

	/* store the proc, offset, and len information in an array
         * of structures, my_req. Each structure contains the
         * offsets and lengths located in that process's FD,
	 * and the associated count.
	 */
	my_req[proc].offsets[l] = off;
  ADIOI_Assert(fd_len == (int) fd_len);
	my_req[proc].lens[l] = (int) fd_len;
	my_req[proc].count++;

	while (rem_len > 0) {
	    off += fd_len;
	    fd_len = rem_len;
	    proc = ADIOI_BG_Calc_aggregator(fd, off, min_st_offset, &fd_len,
					 fd_size, fd_start, fd_end);

	    if (buf_idx[proc] == -1)
      {
        ADIOI_Assert(curr_idx == (int) curr_idx);
        buf_idx[proc] = (int) curr_idx;
      }

	    l = my_req[proc].count;
	    curr_idx += fd_len;
	    rem_len -= fd_len;

	    my_req[proc].offsets[l] = off;
      ADIOI_Assert(fd_len == (int) fd_len);
	    my_req[proc].lens[l] = (int) fd_len;
	    my_req[proc].count++;
	}
    }



#ifdef AGG_DEBUG
    for (i=0; i<nprocs; i++) {
	if (count_my_req_per_proc[i] > 0) {
	    DBG_FPRINTF(stderr, "data needed from %d (count = %d):\n", i,
		    my_req[i].count);
	    for (l=0; l < my_req[i].count; l++) {
		DBG_FPRINTF(stderr, "   off[%d] = %lld, len[%d] = %d\n", l,
			my_req[i].offsets[l], l, my_req[i].lens[l]);
	    }
	}
	DBG_FPRINTF(stderr, "buf_idx[%d] = 0x%x\n", i, buf_idx[i]);
    }
#endif

    *count_my_req_procs_ptr = count_my_req_procs;
    *buf_idx_ptr = buf_idx;
#ifdef AGGREGATION_PROFILE
    MPE_Log_event (5025, 0, NULL);
#endif
    TRACE_ERR("Leaving ADIOI_BG_Calc_my_req\n");
}

/*
 * ADIOI_Calc_others_req (copied to bg and switched to all to all for performance)
 *
 * param[in]  count_my_req_procs        Number of processes whose file domain my
 *                                        request touches.
 * param[in]  count_my_req_per_proc     count_my_req_per_proc[i] gives the no. of
 *                                        contig. requests of this process in
 *                                        process i's file domain.
 * param[in]  my_req                    A structure defining my request
 * param[in]  nprocs                    Number of nodes in the block
 * param[in]  myrank                    Rank of this node
 * param[out] count_others_req_proc_ptr Number of processes whose requests lie in
 *                                        my process's file domain (including my
 *                                        process itself)
 * param[out] others_req_ptr            Array of other process' requests that lie
 *                                        in my process's file domain
 */
void ADIOI_BG_Calc_others_req(ADIO_File fd, int count_my_req_procs,
				int *count_my_req_per_proc,
				ADIOI_Access *my_req,
				int nprocs, int myrank,
				int *count_others_req_procs_ptr,
				ADIOI_Access **others_req_ptr)
{
    TRACE_ERR("Entering ADIOI_BG_Calc_others_req\n");
/* determine what requests of other processes lie in this process's
   file domain */

/* count_others_req_procs = number of processes whose requests lie in
   this process's file domain (including this process itself)
   count_others_req_per_proc[i] indicates how many separate contiguous
   requests of proc. i lie in this process's file domain. */

    int *count_others_req_per_proc, count_others_req_procs;
    int i;
    ADIOI_Access *others_req;

    /* Parameters for MPI_Alltoallv */
    int *scounts, *sdispls, *rcounts, *rdispls;

    /* Parameters for MPI_Alltoallv.  These are the buffers, which
     * are later computed to be the lowest address of all buffers
     * to be sent/received for offsets and lengths.  Initialize to
     * the highest possible address which is the current minimum.
     */
    void *sendBufForOffsets=(void*)0xFFFFFFFFFFFFFFFF,
	 *sendBufForLens   =(void*)0xFFFFFFFFFFFFFFFF,
	 *recvBufForOffsets=(void*)0xFFFFFFFFFFFFFFFF,
	 *recvBufForLens   =(void*)0xFFFFFFFFFFFFFFFF;

/* first find out how much to send/recv and from/to whom */
#ifdef AGGREGATION_PROFILE
    MPE_Log_event (5026, 0, NULL);
#endif
    /* Send 1 int to each process.  count_my_req_per_proc[i] is the number of
     * requests that my process will do to the file domain owned by process[i].
     * Receive 1 int from each process.  count_others_req_per_proc[i] is the number of
     * requests that process[i] will do to the file domain owned by my process.
     */
    count_others_req_per_proc = (int *) ADIOI_Malloc(nprocs*sizeof(int));
/*     cora2a1=timebase(); */
for(i=0;i<nprocs;i++)
    MPI_Alltoall(count_my_req_per_proc, 1, MPI_INT,
		 count_others_req_per_proc, 1, MPI_INT, fd->comm);

/*     total_cora2a+=timebase()-cora2a1; */

    /* Allocate storage for an array of other nodes' accesses of our
     * node's file domain.  Also allocate storage for the alltoallv
     * parameters.
     */
    *others_req_ptr = (ADIOI_Access *)
	ADIOI_Malloc(nprocs*sizeof(ADIOI_Access));
    others_req = *others_req_ptr;

    scounts = ADIOI_Malloc(nprocs*sizeof(int));
    sdispls = ADIOI_Malloc(nprocs*sizeof(int));
    rcounts = ADIOI_Malloc(nprocs*sizeof(int));
    rdispls = ADIOI_Malloc(nprocs*sizeof(int));

    /* If process[i] has any requests in my file domain,
     *   initialize an ADIOI_Access structure that will describe each request
     *   from process[i].  The offsets, lengths, and buffer pointers still need
     *   to be obtained to complete the setting of this structure.
     */
    count_others_req_procs = 0;
    for (i=0; i<nprocs; i++) {
	if (count_others_req_per_proc[i])
  {
	    others_req[i].count = count_others_req_per_proc[i];

	    others_req[i].offsets = (ADIO_Offset *)
		ADIOI_Malloc(count_others_req_per_proc[i]*sizeof(ADIO_Offset));
	    others_req[i].lens = (int *)
		ADIOI_Malloc(count_others_req_per_proc[i]*sizeof(int));

	    if ( (MPIR_Upint)others_req[i].offsets < (MPIR_Upint)recvBufForOffsets )
		recvBufForOffsets = others_req[i].offsets;
	    if ( (MPIR_Upint)others_req[i].lens < (MPIR_Upint)recvBufForLens )
		recvBufForLens = others_req[i].lens;

	    others_req[i].mem_ptrs = (MPI_Aint *)
		ADIOI_Malloc(count_others_req_per_proc[i]*sizeof(MPI_Aint));

	    count_others_req_procs++;
	}
	else
	{
	    others_req[i].count = 0;
	    others_req[i].offsets = NULL;
	    others_req[i].lens    = NULL;
	}
    }
    /* If no recv buffer was allocated in the loop above, make it NULL */
    if ( recvBufForOffsets == (void*)0xFFFFFFFFFFFFFFFF) recvBufForOffsets = NULL;
    if ( recvBufForLens    == (void*)0xFFFFFFFFFFFFFFFF) recvBufForLens    = NULL;

    /* Now send the calculated offsets and lengths to respective processes */

    /************************/
    /* Exchange the offsets */
    /************************/

    /* Determine the lowest sendBufForOffsets/Lens */
    for (i=0; i<nprocs; i++)
    {
	if ( (my_req[i].count) &&
	     ((MPIR_Upint)my_req[i].offsets <= (MPIR_Upint)sendBufForOffsets) )
       {
	  sendBufForOffsets = my_req[i].offsets;
    }

	if ( (my_req[i].count) &&
	     ((MPIR_Upint)my_req[i].lens <= (MPIR_Upint)sendBufForLens) )
       {
	    sendBufForLens = my_req[i].lens;
      }
    }

    /* If no send buffer was found in the loop above, make it NULL */
    if ( sendBufForOffsets == (void*)0xFFFFFFFFFFFFFFFF) sendBufForOffsets = NULL;
    if ( sendBufForLens    == (void*)0xFFFFFFFFFFFFFFFF) sendBufForLens    = NULL;

    /* Calculate the displacements from the sendBufForOffsets/Lens */
    MPI_Barrier(fd->comm);
    for (i=0; i<nprocs; i++)
    {
	/* Send these offsets to process i.*/
	scounts[i] = count_my_req_per_proc[i];
	if ( scounts[i] == 0 )
	    sdispls[i] = 0;
	else
	  sdispls[i] =  (int)
	                ( ( (MPIR_Upint)my_req[i].offsets -
			   (MPIR_Upint)sendBufForOffsets ) /
			  (MPIR_Upint)sizeof(ADIO_Offset) );

	/* Receive these offsets from process i.*/
	rcounts[i] = count_others_req_per_proc[i];
	if ( rcounts[i] == 0 )
	    rdispls[i] = 0;
	else
	    rdispls[i] = (int)
	                 ( ( (MPIR_Upint)others_req[i].offsets -
			     (MPIR_Upint)recvBufForOffsets ) /
			   (MPIR_Upint)sizeof(ADIO_Offset) );
    }

    /* Exchange the offsets */
    MPI_Alltoallv(sendBufForOffsets,
		  scounts, sdispls, ADIO_OFFSET,
		  recvBufForOffsets,
		  rcounts, rdispls, ADIO_OFFSET,
		  fd->comm);

    /************************/
    /* Exchange the lengths */
    /************************/

    for (i=0; i<nprocs; i++)
    {
	/* Send these lengths to process i.*/
	scounts[i] = count_my_req_per_proc[i];
	if ( scounts[i] == 0 )
	    sdispls[i] = 0;
	else
	  sdispls[i] = (int)
	               ( ( (MPIR_Upint)my_req[i].lens -
			   (MPIR_Upint)sendBufForLens ) /
			 (MPIR_Upint) sizeof(int) );

	/* Receive these offsets from process i. */
	rcounts[i] = count_others_req_per_proc[i];
	if ( rcounts[i] == 0 )
	    rdispls[i] = 0;
	else
	    rdispls[i] = (int)
	                 ( ( (MPIR_Upint)others_req[i].lens -
			     (MPIR_Upint)recvBufForLens ) /
			   (MPIR_Upint) sizeof(int) );
    }

    /* Exchange the lengths */
    MPI_Alltoallv(sendBufForLens,
		  scounts, sdispls, MPI_INT,
		  recvBufForLens,
		  rcounts, rdispls, MPI_INT,
		  fd->comm);

    /* Clean up */
    ADIOI_Free(count_others_req_per_proc);
    ADIOI_Free (scounts);
    ADIOI_Free (sdispls);
    ADIOI_Free (rcounts);
    ADIOI_Free (rdispls);

    *count_others_req_procs_ptr = count_others_req_procs;
#ifdef AGGREGATION_PROFILE
    MPE_Log_event (5027, 0, NULL);
#endif
    TRACE_ERR("Leaving ADIOI_BG_Calc_others_req\n");
}
