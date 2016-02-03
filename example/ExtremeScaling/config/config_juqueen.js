P4EST_DIR=/homec/hbn26/hbn264/source/p4est

$P4EST_DIR/configure --host=powerpc64-bgq-linux CFLAGS="-I/bgsys/local/szip/include -I/bgsys/local/zlib/include -shared-libgcc -I/bgsys/drivers/ppcfloor/comm/include -I/bgsys/drivers/ppcfloor/arch/include -I/bgsys/drivers/ppcfloor/gnu-linux/powerpc-bgp-linux/sys-include -O3" CC=mpicc CXX=mpig++ F77=mpigfortran FC=mpigfortran F90=mpigfortran --enable-mpi --disable-shared --without-blas LDFLAGS="-L/bgsys/local/szip/lib -L/bgsys/local/zlib/lib" LIBS="-ldl -lsz"
