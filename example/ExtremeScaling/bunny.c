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

#include <p8est_connectivity.h>
#include <p8est_geometry.h>
#include <p8est_tets_hexes.h>
#include <p8est_vtk.h>

typedef struct
{
    p4est_gloidx_t xm,xM,ym,yM,zm,zM;
}box_t;

static int
bunny_refine (p8est_t * p8est, p4est_topidx_t which_tree,
              p8est_quadrant_t * quadrant)
{
    box_t * box;

    box = (box_t *) p8est->user_pointer;
}

int
main (int argc, char **argv)
{
  int                 mpiret, retval;
  int                 mpirank;
  const char         *argbasename;
  char                afilename[BUFSIZ];
  p4est_topidx_t      tnum_flips;
  p8est_tets_t       *ptg;
  p8est_connectivity_t *connectivity;
  p8est_t            *p8est;
  sc_MPI_Comm         mpicomm;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpicomm = sc_MPI_COMM_WORLD;
  mpiret = sc_MPI_Comm_rank (mpicomm, &mpirank);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  if (argc != 2) {
    SC_GLOBAL_LERRORF ("Usage: %s <tetgen file base name>\n", argv[0]);
    sc_abort ();
  }
  argbasename = argv[1];

  /* read tetgen nodes and tetrahedra from files */
  ptg = p8est_tets_read (argbasename);
  SC_CHECK_ABORTF (ptg != NULL, "Failed to read tetgen %s", argbasename);
  P4EST_GLOBAL_STATISTICSF ("Read %d nodes and %d tets %s attributes\n",
                            (int) ptg->nodes->elem_count / 3,
                            (int) ptg->tets->elem_count / 4,
                            ptg->tet_attributes != NULL ? "with" : "without");

  /* flip orientation to right-handed */
  tnum_flips = p8est_tets_make_righthanded (ptg);
  P4EST_GLOBAL_STATISTICSF ("Performed %ld orientation flip(s)\n",
                            (long) tnum_flips);

  /* create a connectivity from the tet mesh and save it */
  if (mpirank == 0) {
      connectivity = p8est_connectivity_new_tets (ptg);
  }
  else {
      connectivity = NULL;
  }
  //connectivity = p8est_connectivity_bcast (connectivity, 0, sc_MPI_COMM_WORLD);

  P4EST_GLOBAL_LDEBUGF ("Created and broadcasted %s\n", "conn");

  p8est_connectivity_complete (connectivity);

  P4EST_GLOBAL_LDEBUGF ("Connectivity has %ld edges and %ld corners\n",
                        (long) connectivity->num_edges,
                        (long) connectivity->num_corners);

  /*
  if (mpirank == 0) {
    snprintf (afilename, BUFSIZ, "%s", "read_tetgen.p8c");
    retval = p8est_connectivity_save (afilename, connectivity);
    SC_CHECK_ABORT (retval == 0, "Failed connectivity_save");
  }
  */

  /* create a forest and visualize */
  p8est = p8est_new (mpicomm, connectivity, 0, NULL, NULL);
  snprintf (afilename, BUFSIZ, "%s", "read_tetgen");
  p8est_vtk_write_file (p8est, NULL, afilename);

  /* clean up */
  p8est_destroy (p8est);
  p8est_connectivity_destroy (connectivity);
  p8est_tets_destroy (ptg);

  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
