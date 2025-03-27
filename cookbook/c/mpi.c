#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <charm/charm.h>


int main(int argc, char *argv[])
{
    /* Initialize MPI execution environment */
    /* --------------------------------------------------------------------- */
    MPI_Init(&argc, &argv);


    /* Note that the previous "MPI_Init" call does not allow to combined MPI
     * and OpenMP parallelization techniques.  For best performance, however,
     * it is strongly advised to indeed combine MPI and OpenMP.
     *
     * Briefly, the recommended strategy is to use one MPI process per
     * shared-memory system.  This will minimize the data transfer between the
     * computing nodes, which, in case of high-degree spherical harmonics
     * transforms, is usually the bottleneck.  Then, within each shared-memory
     * system, all OpenMP threads can be used.
     *
     * For best performance, we recommend to use the MPI's
     * "MPI_THREAD_FUNNELED" parallelization mode (only the thread that called
     * "MPI_Init_thread" will make MPI calls).  Remember that to combine MPI
     * and OpenMP, the MPI execution environment has to be initialized like
     * this:
     *
     * int provided;
     * MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
     * if (provided != MPI_THREAD_FUNNELED)
     * {
     *     fprintf(stderr, "MPI didn't provide MPI_THREAD_FUNNELED\n");
     *     exit(CHARM_FAILURE);
     * }
     *
     * */
    /* --------------------------------------------------------------------- */


    /* Identify who I am and the number of MPI processes */
    /* --------------------------------------------------------------------- */
    /* MPI communicator */
    MPI_Comm comm = MPI_COMM_WORLD;


    int rank;  /* This is who I am */
    int size;  /* This is the total number of MPI processes */
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);


    /* We assume this code is launched at 3 MPI processes. */
    if (size != 3)
    {
        fprintf(stderr, "This program must be launched at 3 MPI processes.\n");
        exit(CHARM_FAILURE);
    }
    /* --------------------------------------------------------------------- */


    /* Initialize an error structure */
    /* --------------------------------------------------------------------- */
    /* Note that we do not call "charm_err_init" with distributed computing */
    charm_err *err = charm_mpi_err_init(comm);


    /* If the previous function call failed at any process, "err" will be
     * "NULL" at all processes in "comm". */
    if (err == NULL)
    {
        fprintf(stderr, "Failed to create the \"charm_err\" structure");
        exit(CHARM_FAILURE);
    }
    /* --------------------------------------------------------------------- */


    /* Make up some fake spherical harmonic coefficients up to degree "10" */
    /* --------------------------------------------------------------------- */
    /* Maximum harmonic degree of spherical harmonic coefficients */
    unsigned long nmax = 10;


    /* Each MPI process may hold zero, one or more chunks of spherical harmonic
     * coefficients.  A chunk of spherical harmonic coefficients denotes all
     * spherical harmonic coefficients of orders "m1, m1 + 1, ..., m2".
     *
     * In this example, the MPI processes with ranks "0" and "1" will hold two
     * chunks while only a single chunk will be assigned to rank "2".  This is
     * just to demonstrate that the number of chunks may vary across MPI
     * processes.
     *
     * Not shown here, but it is absolutely OK for some processes to store zero
     * chunks, that is, no spherical harmonic coefficients at all.  For
     * instance, if you have "100" MPI processes and you want to create the
     * "charm_shc" structure up to degree "0", you assign one chunk starting
     * and ending at order "0" to one process and zero chunks to all the
     * remaining "99" processes.  */
    size_t local_nchunk;
    if (rank == 2)
        local_nchunk = 1;
    else
        local_nchunk = 2;


    /* To split spherical harmonic coefficients among MPI processes (that is,
     * to define chunks), we need an array of "2 * local_nchunk" elements.  For
     * ranks "0" and "1", we thus need four elements but only two elements for
     * rank "2".  For simplicity, we allocate here, however, an array of four
     * elements for all processes which is fine.
     *
     * In the example below, we hard-coded the chunks for demonstration
     * purposes.  However, this is likely not the strategy you will use in
     * practice, especially if the number of MPI processes is large.  Instead,
     * you'll probably fill "local_order" based on the "local_nchunk", "size"
     * and "rank" variables or so. */
    unsigned long local_order[4];
    if (rank == 0)
    {
        /* The first chunk will consist of all spherical harmonic coefficients
         * of order "0" */
        local_order[0] = 0;
        local_order[1] = 0;


        /* The second chunk will hold all coefficients of orders from "1" to
         * "2" */
        local_order[2] = 1;
        local_order[3] = 2;
    }
    else if (rank == 1)
    {
        /* The first chunk at another MPI process will consist of all
         * coefficients of orders from "3" to "5" */
        local_order[0] = 3;
        local_order[1] = 5;


        /* The second chunk consists of all coefficients of order "6" */
        local_order[2] = 6;
        local_order[3] = 6;
    }
    else if (rank == 2)
    {
        /* This process has one chunk only, starting at order "7" and ending at
         * "10" */
        local_order[0] = 7;
        local_order[1] = 10;
    }


    /* Note that CHarm offers rather flexible distribution of spherical
     * harmonic coefficients among MPI processes.  This has many benefits.  If
     * one of your computing nodes has lower/higher RAM capacity than the other
     * nodes, you can assign lower/higher number of spherical harmonic
     * coefficients to this process.  Or you can distribute all coefficients
     * evenly among all processes to ensure the memory consumption on all nodes
     * will be even and so on and so forth. */


    /* Now allocate the memory for a distributed "charm_shc" structure.  The
     * structure will be distributed among all MPI processed in "comm" based on
     * our "local_nchunk" and "local_order". */
    charm_shc *shcs = charm_mpi_shc_malloc(nmax, 1.0, 1.0, local_nchunk,
                                           local_order, comm, err);


    /* At this point, we have to check that the previous function call was
     * successful.  We can do this either be verifying that "shcs" is not
     * "NULL" and/or we can use the CHarm's error handler.
     *
     * Many things could go wrong when calling "charm_mpi_shc_malloc" and other
     * routines from the "mpi" module.  For instance, there might be not enough
     * memory available at some computing node(s), the partitioning of
     * spherical harmonic coefficients to chunks might be wrong (e.g., by
     * leaving some gaps by not specifying some spherical harmonic order(s) in
     * "local_order", or by introducing overlaps between some chunks), the MPI
     * communicator "comm" might be invalid and so on. You should therefore
     * always call an error handler.  It will try to describe the cause of the
     * error, hopefully allowing you to easily fix it. */
    charm_err_handler(err, 1);


    /* Now fill "shcs" with some fake spherical harmonic coefficients.
     * Usually, you will read these coefficients from files on your hard drive.
     * But for simplicity, we use some artificial coefficients in this example.
     * */
    for (size_t j = 0; j < shcs->local_nchunk; j++)
    {
        for (unsigned long m = shcs->local_order[2 * j];
             m <= shcs->local_order[2 * j + 1]; m++)
        {
            for (unsigned long n = m; n <= shcs->nmax; n++)
            {
                shcs->c[m][n - m] = (double)(m + n);
                shcs->s[m][n - m] = (double)(m * n);
            }
        }
    }


    /* Important: Never access coefficients that are not owned by the process
     * that is trying to access the coefficients.
     *
     * For instance, all coefficients of order "8" are now owned by (stored at)
     * MPI process with rank "2".  This
     *
     *      "shcs->c[8][8 - 8]"
     *
     * will thus fail badly at processes with ranks other than "2". */
    /* --------------------------------------------------------------------- */


    /* Now let's prepare a Gauss--Legendre grid that will be distributed among
     * the MPI processes */
    /* --------------------------------------------------------------------- */
    /* Get the number of latitudes and longitudes of the Gauss--Legendre grid
     * that correspond to maximum harmonic degree "nmax" */
    size_t nlat, nlon;
    charm_crd_point_gl_shape(nmax, &nlat, &nlon);


    /* For our "nmax = 10", "nlat" will be "11".  We will thus have the
     * following latitudes in the Gauss--Legendre grid: "lat0, lat1, lat2, ...,
     * lat10". */


    /* Similarly as we specified chunks of spherical harmonic coefficients, now
     * we have to split the Gauss--Legendre grid into latitudinal chunks.
     * Naturally, the sum of "local_nlat" across all processes must match
     * "nlat".  This is automatically checked when calling
     * "charm_mpi_crd_point_gl" below.
     *
     * Again, assigning no latitudinal chunks to some MPI process is fairly OK.
     * Below, this would translate into assigning "0" to "local_nlat" at one or
     * more processes.  This would of course necessitate to increase
     * "local_nlat" at other processes. */
    size_t local_nlat, local_0_start;
    if (rank == 0)
    {
        /* This MPI process will hold 2 latitudes, "lat1" and "lat9". */
        local_nlat    = 2;
        local_0_start = 1;
    }
    else if (rank == 1)
    {
        /* This process also holds two latitudes, "lat0" and "lat10".  We
         * intentionally put "lat0" to the process of rank "1" to demonstrate
         * that CHarm does not care about which chunk is stored by which MPI
         * process. */
        local_nlat    = 2;
        local_0_start = 0;
    }
    else if (rank == 2)
    {
        /* This process holds seven latitudes, "lat2, lat3, ..., lat8" */
        local_nlat    = 7;
        local_0_start = 2;
    }
    /* Note that "lat0 = -lat10", "lat1 = -lat9", ..., that is Gauss--Legendre
     * grids are symmetric with respect to the equator.  In case of symmetric
     * grids in general, each process must hold the latitudes at both
     * hemispheres.  The only exception is that MPI process, which holds the
     * equator (if there is any).  In that case, "local_nlat" must be odd.
     * This happens with Gauss--Legendre grids and even "nmax".  See the
     * documentation to "charm_mpi_crd_point_gl". */


    charm_point *glg = charm_mpi_crd_point_gl(nmax, shcs->r, local_nlat,
                                              local_0_start, comm, err);
    charm_err_handler(err, 1);
    /* --------------------------------------------------------------------- */


    /* Let's print the latitudes, spherical radii and integration weights in
     * "glg".  The longitudes are not printed as they are the same at all
     * processes.  The data will be printed in the rank order. */
    /* --------------------------------------------------------------------- */
    for (int r = 0; r < size; r++)
    {
        if (rank == r)
        {
            for (size_t i = 0; i < glg->local_nlat; i++)
                printf("rank %d: glg->lat[%zu] = %23.16e, "
                       "glg->r[%zu] = %23.16e, glg->w[%zu] = %23.16e\n",
                       rank, i, glg->lat[i], i, glg->r[i], i, glg->w[i]);
        }


        MPI_Barrier(comm);
    }
    /* --------------------------------------------------------------------- */


    /* At this point, we have spherical harmonic coefficients in "shcs" as well
     * as computing points in "glg" distributed across MPI processes, so we can
     * do some distributed spherical harmonic synthesis and analysis. */


    /* Spherical harmonic synthesis */
    /* --------------------------------------------------------------------- */
    /* At this point, each process holds "glg->local_npoint" points.  The sum
     * of "glg->local_npoint" across all MPI processes thus yields the total
     * number of grid points "glg->npoint".  Let's now allocate memory for the
     * local signal "f" that is to be synthesized by each MPI process */
    double *f = (double *)malloc(glg->local_npoint * sizeof(double));
    if (f == NULL)
    {
        fprintf(stderr, "malloc failure.\n");
        exit(CHARM_FAILURE);
    }


    /* Do the synthesis */
    charm_shs_point(glg, shcs, shcs->nmax, f, err);
    charm_err_handler(err, 1);


    /* Note that "f" is now distributed among MPI processes.  The distribution
     * is given by "glg". */
    /* --------------------------------------------------------------------- */


    /* Spherical harmonic analysis */
    /* --------------------------------------------------------------------- */
    /* Now we want to compute by spherical harmonic analysis coefficients from
     * "f", which is distributed.  First, we need to allocate a new distributed
     * "charm_shc" structure.  For simplicity, the distribution will be the
     * same as with "shcs" but this is not required. */
    charm_shc *shcs_out = charm_mpi_shc_malloc(nmax, 1.0, 1.0, local_nchunk,
                                               local_order, comm, err);
    charm_err_handler(err, 1);


    /* Now the analysis.  Note that the number of locally stored points in
     * "glg" ("glg->local_npoint") must match the number of elements in "f" and
     * that the elements of "f" must be ordered properly. */
    charm_sha_point(glg, f, shcs_out->nmax, shcs_out, err);
    charm_err_handler(err, 1);


    /* "shcs" and "shcs_out" are the same, as can easily be verified.  Their
     * differences should not exceed the order of, say, "10**-13" or so.  We
     * print the outputs first in order according to the rank of the process
     * and then in order given by "local_order". */
    for (int r = 0; r < size; r++)
    {
        if (rank == r)
        {
            for (size_t j = 0; j < shcs->local_nchunk; j++)
            {
                for (unsigned long m = shcs->local_order[2 * j];
                     m <= shcs->local_order[2 * j + 1]; m++)
                {
                    for (unsigned long n = m; n <= shcs->nmax; n++)
                    {
                        printf("rank %d: Cref_{%2lu, %2lu} = %0.16e, "
                               "Cnew_{%2lu, %2lu} = %0.16e, diff = %0.16e\n",
                               rank, n, m, shcs->c[m][n - m], n, m,
                               shcs_out->c[m][n - m],
                               shcs->c[m][n - m] - shcs_out->c[m][n - m]);
                        printf("rank %d: Sref_{%2lu, %2lu} = %0.16e, "
                               "Snew_{%2lu, %2lu} = %0.16e, diff = %0.16e\n",
                               rank, n, m, shcs->s[m][n - m], n, m,
                               shcs_out->s[m][n - m],
                               shcs->s[m][n - m] - shcs_out->s[m][n - m]);
                    }
                }
            }
        }


        MPI_Barrier(comm);
    }
    /* --------------------------------------------------------------------- */


    /* IMPORTANT NOTE */
    /* --------------------------------------------------------------------- */
    /* CHarm has some global parameters that should be tuned to get the best
     * performance with MPI.  For instance, there is the
     * "charm_glob_sha_block_lat_multiplier" variable, which has enormous
     * impact on the performance of "charm_sha_point" in real-world high-degree
     * applications.
     *
     * For demonstration purposes, let's change its value now and repeat the
     * analysis.  This has no effect on the results of course.
     *
     * See the documentation of the "glob" module for this and other associated
     * variables. */
    charm_glob_sha_block_lat_multiplier = 100;
    charm_sha_point(glg, f, shcs_out->nmax, shcs_out, err);
    charm_err_handler(err, 1);
    /* --------------------------------------------------------------------- */


    /* Free the memory */
    /* --------------------------------------------------------------------- */
    free(f);
    charm_shc_free(shcs);
    charm_shc_free(shcs_out);
    charm_crd_point_free(glg);
    charm_err_free(err);
    /* --------------------------------------------------------------------- */


    /* --------------------------------------------------------------------- */
    printf("Great, all done by process %d of %d!\n", rank, size);


    MPI_Finalize();
    /* --------------------------------------------------------------------- */


    exit(CHARM_SUCCESS);
    /* ===================================================================== */
}
