/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2010 The University of Texas System
  Written by Carsten Burstedde, Lucas C. Wilcox, and Tobin Isaac

  p4est is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  p4est is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with p4est; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#ifndef P4_TO_P8
#include <p4est_extended.h>
#include <p4est_ghost.h>
#else
#include <p8est_extended.h>
#include <p8est_geometry.h>
#include <p8est_ghost.h>
#include <p8est_tets_hexes.h>
#endif

#include <sc_flops.h>
#include <sc_statistics.h>
#include <sc_options.h>

typedef enum exchange_timers
{
  EXCHANGE_ZERO = 0,
  EXCHANGE_READTETS = EXCHANGE_ZERO,
  EXCHANGE_HANDED,
  EXCHANGE_FROMTETS,
  EXCHANGE_COMPLETE,
  EXCHANGE_BCAST,
  EXCHANGE_P4EST,
  EXCHANGE_REFINE,
  EXCHANGE_PARTITION,
  EXCHANGE_GHOST,
  EXCHANGE_FILL,
  EXCHANGE_EXCHANGE,
  EXCHANGE_VERIFY,
  EXCHANGE_TIMERS
}
exchange_timers_t;

static const char  *timer_names[EXCHANGE_TIMERS] = {
  "Readtets", "Handed", "Fromtets", "Complete",
  "Bcast", "P4est", "Refine", "Partition",
  "Ghost", "Fill", "Exchange", "Verify"
};

typedef struct
{
  double              xm, xM, ym, yM, zm, zM;
}
box_t;

#ifdef P4_TO_P8

/** Get the coordinates of the midpoint of a quadrant.
 *
 * \param [in]  p4est      the forest
 * \param [in]  which_tree the tree in the forest containing \a q
 * \param [in]  q          the quadrant
 * \param [out] xyz        the coordinates of the midpoint of \a q
 */
static void
bunny_get_midpoint (p4est_t * p4est, p4est_topidx_t which_tree,
                    p4est_quadrant_t * q, double xyz[3])
{
  p4est_qcoord_t      half_length = P4EST_QUADRANT_LEN (q->level) / 2;

  p4est_qcoord_to_vertex (p4est->connectivity, which_tree,
                          q->x + half_length, q->y + half_length,
                          q->z + half_length, xyz);
}

/* Refine if we lie in a cone defined by a bounding box */
static int
bunny_refine (p4est_t * p4est, p4est_topidx_t which_tree,
              p4est_quadrant_t * quadrant)
{
  box_t              *box;
  double              coords[3];
  double              R, r;

  box = (box_t *) p4est->user_pointer;

  bunny_get_midpoint (p4est, which_tree, quadrant, coords);
  R = (box->xM - box->xm) / 4.;
  r = (coords[1] - box->ym) / (box->yM - box->ym) * R;
  if (pow ((coords[0] - (box->xM + box->xm) / 2), 2)
      + pow ((coords[2] - (box->zM + box->zm) / 2), 2) <= r * r) {
    return 1;
  }
  else
    return 0;
}

#endif /* P4_TO_P8 */

static const p4est_locidx_t encode = 0xdead13af;

static void
fill_ghost_data (p4est_t * p4est, p4est_ghost_t * ghost)
{
#ifdef P4EST_ENABLE_DEBUG
  size_t              data_size;
#endif
  size_t              zz;
  p4est_topidx_t      tree_id;
  p4est_locidx_t      local_id, intree_pos;
  p4est_quadrant_t   *mirror, *quadrant;
  p4est_tree_t       *tree;

#ifdef P4EST_ENABLE_DEBUG
  data_size = p4est->data_size;
  P4EST_ASSERT (data_size == sizeof (p4est_locidx_t));
#endif

  for (zz = 0; zz < ghost->mirrors.elem_count; ++zz) {
    /* get tree and local quadrant number from mirror data */
    mirror = p4est_quadrant_array_index (&ghost->mirrors, zz);
    tree_id = mirror->p.piggy3.which_tree;
    P4EST_ASSERT (p4est->first_local_tree <= tree_id &&
                  tree_id <= p4est->last_local_tree);
    local_id = mirror->p.piggy3.local_num;
    P4EST_ASSERT (0 <= local_id && local_id < p4est->local_num_quadrants);

    /* figure out the real existing quadrant pointed at by the mirror */
    tree = p4est_tree_array_index (p4est->trees, tree_id);
    P4EST_ASSERT (tree->quadrants_offset <= local_id);
    intree_pos = local_id - tree->quadrants_offset;
    P4EST_ASSERT (intree_pos < (p4est_locidx_t) tree->quadrants.elem_count);
    quadrant = p4est_quadrant_array_index (&tree->quadrants, intree_pos);

    /* assign this number into the data field and encode it */
    *(p4est_locidx_t *) quadrant->p.user_data = local_id ^ encode;
  }
}

static void
verify_ghost_data (p4est_t * p4est, p4est_ghost_t * ghost,
                   p4est_locidx_t * received)
{
#ifdef P4EST_ENABLE_DEBUG
  size_t              data_size;
#endif
  size_t              zz;
#if 0
  p4est_topidx_t      tree_id;
#endif
  p4est_locidx_t      remote_id;
  p4est_quadrant_t   *remote;

#ifdef P4EST_ENABLE_DEBUG
  data_size = p4est->data_size;
  P4EST_ASSERT (data_size == sizeof (p4est_locidx_t));
#endif

  for (zz = 0; zz < ghost->ghosts.elem_count; ++zz) {
    /* get tree and remote quadrant number from ghost data */
    remote = p4est_quadrant_array_index (&ghost->ghosts, zz);
#if 0
    tree_id = remote->p.piggy3.which_tree;
    P4EST_ASSERT (0 <= tree_id && tree_id < p4est->connectivity->num_trees);
#endif
    remote_id = remote->p.piggy3.local_num ^ encode;

    SC_CHECK_ABORT (remote_id == received[zz], "Ghost data mismatch");
  }
}

int
main (int argc, char **argv)
{
#if 0
  int                 retval;
  char                afilename[BUFSIZ];
#endif
  int                 mpiret;
  int                 mpirank;
  int                 irun, runs;
  int                 extim;
  size_t              data_size;
  double              snaptime[EXCHANGE_TIMERS];
  p4est_locidx_t     *received;
  sc_MPI_Comm         mpicomm;
  sc_flopinfo_t       fi, snapshot;
  sc_statinfo_t       stats[EXCHANGE_TIMERS];
  p4est_connectivity_t *connectivity;
  p4est_t            *p4est;
  p4est_ghost_t      *ghost;
#ifdef P4_TO_P8
  const char         *argbasename;
  p4est_topidx_t      tnum_flips;
  p8est_tets_t       *ptg;
#endif
  box_t               Box_ex1;

/*
 * In 2D we run moebius.
 * In 3D we run rotcubes and a tetgen input file.
 */

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_STATISTICS);

  mpicomm = sc_MPI_COMM_WORLD;
  mpiret = sc_MPI_Comm_rank (mpicomm, &mpirank);
  SC_CHECK_MPI (mpiret);

#ifndef P4_TO_P8
  runs = 1;
#else
  runs = 2;
#endif

  for (irun = 0; irun < runs; ++irun) {
    P4EST_GLOBAL_ESSENTIALF ("Run %s exchange %d\n", P4EST_STRING, irun);
#ifdef P4_TO_P8
    argbasename = NULL;
    ptg = NULL;
    tnum_flips = 0;
    if (irun == 1) {
      if (argc != 2) {
        SC_GLOBAL_LERRORF ("Usage: %s <tetgen file base name>\n", argv[0]);
        sc_abort ();
      }
      argbasename = argv[1];
    }
#endif

    data_size = sizeof (p4est_locidx_t);

    Box_ex1.xm = -6;
    Box_ex1.ym = -6;
    Box_ex1.zm = -6;
    Box_ex1.xM = 7;
    Box_ex1.yM = 7;
    Box_ex1.zM = 7;

    /* prepare timers */
    for (extim = EXCHANGE_ZERO; extim < EXCHANGE_TIMERS; ++extim) {
      snaptime[extim] = 0.;
    }
    sc_flops_start (&fi);

    connectivity = NULL;
#ifndef P4_TO_P8
    if (mpirank == 0) {
      connectivity = p4est_connectivity_new_moebius ();
    }
#else
    if (irun == 0) {
      if (mpirank == 0) {
        connectivity = p8est_connectivity_new_rotcubes ();
      }
    }
    else {
      P4EST_ASSERT (irun == 1);
      /* so only here we read a mesh file and only on one processor. */

      /* read tetgen nodes and tetrahedra from files */
      sc_flops_snap (&fi, &snapshot);
      if (mpirank == 0) {
        P4EST_ASSERT (argbasename != NULL);
        ptg = p8est_tets_read (argbasename);
        SC_CHECK_ABORTF (ptg != NULL, "Failed to read tetgen %s",
                         argbasename);

        P4EST_GLOBAL_STATISTICSF ("Read %d nodes and %d tets %s attributes\n",
                                  (int) ptg->nodes->elem_count / 3,
                                  (int) ptg->tets->elem_count / 4,
                                  ptg->tet_attributes !=
                                  NULL ? "with" : "without");
      }
      sc_flops_shot (&fi, &snapshot);
      snaptime[EXCHANGE_READTETS] = snapshot.iwtime;

      /* flip orientation to right-handed */
      sc_flops_snap (&fi, &snapshot);
      if (mpirank == 0) {
        tnum_flips = p8est_tets_make_righthanded (ptg);
        P4EST_GLOBAL_STATISTICSF ("Performed %ld orientation flip(s)\n",
                                  (long) tnum_flips);
      }
      sc_flops_shot (&fi, &snapshot);
      snaptime[EXCHANGE_HANDED] = snapshot.iwtime;

      /* create a connectivity from the tet mesh */
      sc_flops_snap (&fi, &snapshot);
      if (mpirank == 0) {
        connectivity = p8est_connectivity_new_tets (ptg);
      }
      sc_flops_shot (&fi, &snapshot);
      snaptime[EXCHANGE_FROMTETS] = snapshot.iwtime;

      /* complete connectivity information that is missing in file */
      sc_flops_snap (&fi, &snapshot);
      if (mpirank == 0) {
        p8est_connectivity_complete (connectivity);

        P4EST_GLOBAL_STATISTICSF
          ("Connectivity has %ld edges and %ld corners\n",
           (long) connectivity->num_edges, (long) connectivity->num_corners);
      }
      sc_flops_shot (&fi, &snapshot);
      snaptime[EXCHANGE_COMPLETE] = snapshot.iwtime;

      /* TODO: save connectivity to p4est binary format */
#if 0
      if (mpirank == 0) {
        snprintf (afilename, BUFSIZ, "%s", "read_tetgen.p8c");
        retval = p8est_connectivity_save (afilename, connectivity);
        SC_CHECK_ABORT (retval == 0, "Failed connectivity_save");
      }
#endif /* 0 */
    }
#endif /* P4_TO_P8 */
    P4EST_ASSERT ((mpirank == 0) == (connectivity != NULL));

    /* this code is the same for all examples */
    sc_flops_snap (&fi, &snapshot);
    connectivity = p4est_connectivity_bcast (connectivity, 0, mpicomm);
    P4EST_GLOBAL_LDEBUGF ("Created and broadcasted %s\n", "conn");
    sc_flops_shot (&fi, &snapshot);
    snaptime[EXCHANGE_BCAST] = snapshot.iwtime;

    /* create a forest */
    sc_flops_snap (&fi, &snapshot);
    p4est = p4est_new_ext (mpicomm, connectivity, 0, 4, 1,
                           data_size, NULL, &Box_ex1);
    sc_flops_shot (&fi, &snapshot);
    snaptime[EXCHANGE_P4EST] = snapshot.iwtime;

    /* refine the forest */
    sc_flops_snap (&fi, &snapshot);
#ifdef P4_TO_P8
    p4est_refine (p4est, 0, bunny_refine, NULL);
#endif
    sc_flops_shot (&fi, &snapshot);
    snaptime[EXCHANGE_REFINE] = snapshot.iwtime;

    /* partition the forest */
    sc_flops_snap (&fi, &snapshot);
    p4est_partition (p4est, 1, NULL);
    sc_flops_shot (&fi, &snapshot);
    snaptime[EXCHANGE_PARTITION] = snapshot.iwtime;

    /* create a ghost layer */
    sc_flops_snap (&fi, &snapshot);
    ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FULL);
    sc_flops_shot (&fi, &snapshot);
    snaptime[EXCHANGE_GHOST] = snapshot.iwtime;

    /* fill data to prepare ghost exchange */
    sc_flops_snap (&fi, &snapshot);
    received = P4EST_ALLOC (p4est_locidx_t, ghost->ghosts.elem_count);
    fill_ghost_data (p4est, ghost);
    sc_flops_shot (&fi, &snapshot);
    snaptime[EXCHANGE_FILL] = snapshot.iwtime;

    /* execute ghost exchange */
    sc_flops_snap (&fi, &snapshot);
    p4est_ghost_exchange_data (p4est, ghost, received);
    sc_flops_shot (&fi, &snapshot);
    snaptime[EXCHANGE_EXCHANGE] = snapshot.iwtime;

    /* verify data after ghost exchange */
    sc_flops_snap (&fi, &snapshot);
    verify_ghost_data (p4est, ghost, received);
    sc_flops_shot (&fi, &snapshot);
    snaptime[EXCHANGE_VERIFY] = snapshot.iwtime;

    /* clean up */
    P4EST_FREE (received);
    p4est_ghost_destroy (ghost);
    p4est_destroy (p4est);
    p4est_connectivity_destroy (connectivity);

#ifdef P4_TO_P8
    if (irun == 1) {
      p8est_tets_destroy (ptg);
    }
#endif

    for (extim = EXCHANGE_ZERO; extim < EXCHANGE_TIMERS; ++extim) {
      sc_stats_set1 (&stats[extim], snaptime[extim], timer_names[extim]);
    }
    sc_stats_compute (mpicomm, EXCHANGE_TIMERS, stats);
    sc_stats_print (p4est_package_id, SC_LP_STATISTICS, EXCHANGE_TIMERS,
                    stats, 1, 1);
  }

  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
