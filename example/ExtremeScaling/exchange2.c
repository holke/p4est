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
#else
#include <p8est_extended.h>
#include <p8est_geometry.h>
#include <p8est_tets_hexes.h>
#endif

#include <sc_flops.h>
#include <sc_statistics.h>
#include <sc_options.h>

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
bunny_get_midpoint (p8est_t * p8est, p4est_topidx_t which_tree,
                    p8est_quadrant_t * q, double xyz[3])
{
  p4est_qcoord_t      half_length = P8EST_QUADRANT_LEN (q->level) / 2;

  p8est_qcoord_to_vertex (p8est->connectivity, which_tree,
                          q->x + half_length, q->y + half_length,
                          q->z + half_length, xyz);
}

/* Refine if we lie in a cylinder defined by a bounding box */
static int
bunny_refine (p8est_t * p8est, p4est_topidx_t which_tree,
              p8est_quadrant_t * quadrant)
{
  box_t              *box;
  double              coords[3];
  double              R, r;

  box = (box_t *) p8est->user_pointer;

  bunny_get_midpoint (p8est, which_tree, quadrant, coords);
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
  double              time_read, time_handed, time_newtets, time_complete;
  p4est_connectivity_t *connectivity;
  p4est_t            *p8est;
  sc_MPI_Comm         mpicomm;
  box_t               Box_ex1;
  sc_flopinfo_t       fi, snapshot;
  sc_statinfo_t       stats[7];
#ifdef P4_TO_P8
  const char         *argbasename;
  p4est_topidx_t      tnum_flips;
  p8est_tets_t       *ptg;
#endif

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

    Box_ex1.xm = -6;
    Box_ex1.ym = -6;
    Box_ex1.zm = -6;
    Box_ex1.xM = 7;
    Box_ex1.yM = 7;
    Box_ex1.zM = 7;

    time_read = time_handed = time_newtets = time_complete = 0.;

    sc_flops_start (&fi);
    sc_flops_snap (&fi, &snapshot);

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

      if (mpirank == 0) {
        /* read tetgen nodes and tetrahedra from files */
        P4EST_ASSERT (argbasename != NULL);
        ptg = p8est_tets_read (argbasename);
        SC_CHECK_ABORTF (ptg != NULL, "Failed to read tetgen %s",
                         argbasename);
      }
      sc_flops_shot (&fi, &snapshot);
      time_read = snapshot.iwtime;
      sc_flops_snap (&fi, &snapshot);

      if (mpirank == 0) {
        P4EST_GLOBAL_STATISTICSF ("Read %d nodes and %d tets %s attributes\n",
                                  (int) ptg->nodes->elem_count / 3,
                                  (int) ptg->tets->elem_count / 4,
                                  ptg->tet_attributes !=
                                  NULL ? "with" : "without");

        /* flip orientation to right-handed */
        tnum_flips = p8est_tets_make_righthanded (ptg);
        P4EST_GLOBAL_STATISTICSF ("Performed %ld orientation flip(s)\n",
                                  (long) tnum_flips);
      }
      sc_flops_shot (&fi, &snapshot);
      time_handed = snapshot.iwtime;
      sc_flops_snap (&fi, &snapshot);

      if (mpirank == 0) {
        /* create a connectivity from the tet mesh */
        connectivity = p8est_connectivity_new_tets (ptg);
      }
      sc_flops_shot (&fi, &snapshot);
      time_newtets = snapshot.iwtime;
      sc_flops_snap (&fi, &snapshot);

      if (mpirank == 0) {
        /* complete connectivity information that is missing in file */
        p8est_connectivity_complete (connectivity);

        P4EST_GLOBAL_LDEBUGF ("Connectivity has %ld edges and %ld corners\n",
                              (long) connectivity->num_edges,
                              (long) connectivity->num_corners);
      }
      sc_flops_shot (&fi, &snapshot);
      time_complete = snapshot.iwtime;
      sc_flops_snap (&fi, &snapshot);

      /* TODO: save connectivity to p4est binary format */
      /*
         if (mpirank == 0) {
         snprintf (afilename, BUFSIZ, "%s", "read_tetgen.p8c");
         retval = p8est_connectivity_save (afilename, connectivity);
         SC_CHECK_ABORT (retval == 0, "Failed connectivity_save");
         }
       */
    }
#endif
    P4EST_ASSERT ((mpirank == 0) == (connectivity != NULL));

    /* this code is the same for all examples */
    connectivity = p4est_connectivity_bcast (connectivity, 0, mpicomm);
    sc_flops_shot (&fi, &snapshot);
    sc_stats_set1 (&stats[2], snapshot.iwtime, "Bcast");
    sc_flops_snap (&fi, &snapshot);

    P4EST_GLOBAL_LDEBUGF ("Created and broadcasted %s\n", "conn");

    /* create a forest and visualize */

    p8est = p4est_new_ext (mpicomm, connectivity, 0, 4, 1, 0, NULL, &Box_ex1);
    sc_flops_shot (&fi, &snapshot);
    sc_stats_set1 (&stats[4], snapshot.iwtime, P4EST_STRING "_new");
    sc_flops_snap (&fi, &snapshot);

#ifdef P4_TO_P8
    p4est_refine (p8est, 0, bunny_refine, NULL);
#endif
    sc_flops_shot (&fi, &snapshot);
    sc_stats_set1 (&stats[5], snapshot.iwtime, "Refine 1 times");
    sc_flops_snap (&fi, &snapshot);

    /* clean up */
    p4est_destroy (p8est);
    p4est_connectivity_destroy (connectivity);

#ifdef P4_TO_P8
    if (irun == 1) {
      p8est_tets_destroy (ptg);
    }
#endif

    sc_stats_set1 (&stats[0], time_read, "Read");
    sc_stats_set1 (&stats[1], time_handed, "Right handed");
    sc_stats_set1 (&stats[3], time_complete, "Connectivity complete");
    sc_stats_set1 (&stats[6], time_newtets, "New connectivity from tets");

    sc_stats_compute (mpicomm, 7, stats);
    sc_stats_print (p4est_package_id, SC_LP_STATISTICS, 7, stats, 1, 1);
  }

  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
