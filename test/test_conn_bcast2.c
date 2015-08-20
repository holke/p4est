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
#include <p4est_connectivity.h>
#else
#include <p8est_connectivity.h>
#endif

static void
test_bcast (p4est_connectivity_t * conn, const char * which, sc_MPI_Comm comm)
{
  int                mpirank, mpisize, mpiret;
  p4est_connectivity_t *conn_recv;

  SC_GLOBAL_INFOF ("Testing standard connectivity %s\n", which);
  SC_CHECK_ABORTF (p4est_connectivity_is_valid (conn),
                   "Invalid connectivity %s before broadcast", which);
  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);
  /* Processor 0 will send the connectivity conn to all other processors.
   * The received conn will be stored in conn_recv and we than check
   * whether conn_recv is valid and the same connectivity as conn */
  if (mpirank == 0) {
    conn_recv = conn;
  }
  else {
    conn_recv = NULL;
  }
  conn_recv = p4est_connectivity_bcast (conn_recv, 0, comm);
  SC_CHECK_ABORTF (p4est_connectivity_is_valid (conn),
                   "Invalid connectivity %s after broadcast", which);
  SC_CHECK_ABORTF (p4est_connectivity_is_equal (conn, conn_recv),
                   "Received connectivity %s does not match the connectivity"
                   "that was sent", which);
  p4est_connectivity_destroy (conn);
  if (mpirank != 0) {
    p4est_connectivity_destroy (conn_recv);
  }
}

int
main (int argc, char *argv[])
{
  int           mpiret;
  sc_MPI_Comm   comm;

  /* initialize MPI and p4est internals */
  comm = sc_MPI_COMM_WORLD;
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  sc_init (comm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

#ifndef P4_TO_P8
  test_bcast (p4est_connectivity_new_unitsquare (), "unitsquare", comm);
  test_bcast (p4est_connectivity_new_periodic (), "periodic", comm);
  test_bcast (p4est_connectivity_new_rotwrap (), "rotwrap", comm);
  test_bcast (p4est_connectivity_new_corner (), "corner", comm);
  test_bcast (p4est_connectivity_new_pillow (), "pillow", comm);
  test_bcast (p4est_connectivity_new_moebius (), "moebius", comm);
  test_bcast (p4est_connectivity_new_star (), "star", comm);
  test_bcast (p4est_connectivity_new_cubed (), "cubed", comm);
  test_bcast (p4est_connectivity_new_disk (), "disk", comm);
  test_bcast (p4est_connectivity_new_brick (3, 2, 0, 0), "brick00", comm);
  test_bcast (p4est_connectivity_new_brick (3, 2, 0, 1), "brick01", comm);
  test_bcast (p4est_connectivity_new_brick (3, 2, 1, 0), "brick10", comm);
  test_bcast (p4est_connectivity_new_brick (3, 2, 1, 1), "brick11", comm);
#else
  test_bcast (p8est_connectivity_new_unitcube (), "unitcube", comm);
  test_bcast (p8est_connectivity_new_periodic (), "periodic", comm);
  test_bcast (p8est_connectivity_new_rotwrap (), "rotwrap", comm);
  test_bcast (p8est_connectivity_new_twocubes (), "twocubes", comm);
  test_bcast (p8est_connectivity_new_twowrap (), "twowrap", comm);
  test_bcast (p8est_connectivity_new_rotcubes (), "rotcubes", comm);
  test_bcast (p8est_connectivity_new_brick (4, 3, 2, 0, 0, 0), "brick000", comm);
  test_bcast (p8est_connectivity_new_brick (4, 3, 2, 0, 0, 1), "brick001", comm);
  test_bcast (p8est_connectivity_new_brick (4, 3, 2, 0, 1, 0), "brick010", comm);
  test_bcast (p8est_connectivity_new_brick (4, 3, 2, 0, 1, 1), "brick011", comm);
  test_bcast (p8est_connectivity_new_brick (4, 3, 2, 1, 0, 0), "brick100", comm);
  test_bcast (p8est_connectivity_new_brick (4, 3, 2, 1, 0, 1), "brick101", comm);
  test_bcast (p8est_connectivity_new_brick (4, 3, 2, 1, 1, 1), "brick111", comm);
#endif
  /* clean up and exit */
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
