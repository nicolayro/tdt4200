#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <unistd.h>
#include <math.h>
#include <mpi.h>
#include <sys/time.h>

#include "../inc/argument_utils.h"

// Convert 'struct timeval' into seconds in double prec. floating point
#define WALLTIME(t) ((double)(t).tv_sec + 1e-6 * (double)(t).tv_usec)

typedef int64_t int_t;
typedef double real_t;

int_t
    N,
    M,
    max_iteration,
    snapshot_frequency;

real_t
    *temp[2] = { NULL, NULL },
    *thermal_diffusivity,
    dx,
    dt;

int 
    size,
    rank,
    len,
    left_neighbor,
    right_neighbor,
    sub_grid_offset;

#define T(i,j)                      temp[0][(i) * (M + 2) + (j)]
#define T_next(i,j)                 temp[1][((i) * (M + 2) + (j))]
#define THERMAL_DIFFUSIVITY(i,j)    thermal_diffusivity[(i) * (M + 2) + (j)]

void time_step ( void );
void boundary_condition( void );
void border_exchange( void );
void domain_init ( void );
void domain_save ( int_t iteration );
void domain_finalize ( void );


void
swap ( real_t** m1, real_t** m2 )
{
    real_t* tmp;
    tmp = *m1;
    *m1 = *m2;
    *m2 = tmp;
}


int
main ( int argc, char **argv )
{
    // TODO 1: Initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &size );

    // TODO 2: Parse arguments in the rank 0 processes
    // and broadcast to other processes
    OPTIONS *options = parse_args( argc, argv );
    if ( !options )
    {
        fprintf( stderr, "Argument parsing failed\n" );
        exit(1);
    }

    // Task 2
    if (rank == 0) {
        N = options->N;
        M = options->M;
        max_iteration = options->max_iteration;
        snapshot_frequency = options->snapshot_frequency;
        for (int i = 1; i < size; i++) {
            MPI_Send(&M, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(&N, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(&max_iteration, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(&snapshot_frequency, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        }
    } else {
        MPI_Recv(&M, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&N, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&max_iteration, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&snapshot_frequency, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    // TODO 3: Allocate space for each process' sub-grids
    // and initialize data for the sub-grids
    domain_init();

    struct timeval t_start, t_end;
    gettimeofday ( &t_start, NULL );

    for ( int_t iteration = 0; iteration <= max_iteration; iteration++ )
    {
        // TODO 7: Communicate border values
        border_exchange();

        // TODO 5: Boundary conditions
        boundary_condition();

        // TODO 4: Time step calculations
        time_step();

        if ( iteration % snapshot_frequency == 0 )
        {
            printf (
                "Iteration %ld of %ld (%.2lf%% complete)\n",
                (long) iteration,
                (long) max_iteration,
                100.0 * (real_t) iteration / (real_t) max_iteration
            );

            // TODO 6 MPI I/O
            domain_save ( iteration );
        }

        swap( &temp[0], &temp[1] );
    }
    gettimeofday ( &t_end, NULL );
    printf ( "Total elapsed time: %lf seconds\n",
            WALLTIME(t_end) - WALLTIME(t_start)
            );

    domain_finalize();

    // TODO 1: Finalize MPI
    MPI_Finalize();

    exit ( EXIT_SUCCESS );
}


void
time_step ( void )
{
    // TODO 4: Time step calculations
    real_t c, t, b, l, r, K, new_value;

    for ( int_t x = 1; x <= N / size; x++ )
    {
        for ( int_t y = 1; y <= M; y++ )
        {
            c = T(x, y);

            t = T(x - 1, y);
            b = T(x + 1, y);
            l = T(x, y - 1);
            r = T(x, y + 1);
            K = THERMAL_DIFFUSIVITY(x, y);

            new_value = c + K * dt * ((l - 2 * c + r) + (b - 2 * c + t));

            T_next(x, y) = new_value;
        }
    }
}


void
boundary_condition ( void )
{
    // TODO 5: Boundary conditions
    for ( int_t x = 1; x <= N / size; x++ )
    {
        T(x, 0) = T(x, 2);
        T(x, M+1) = T(x, M-1);
    }

    for ( int_t y = 1; y <= M; y++ )
    {
        T(0, y) = T(2, y);
        T(N / size + 1, y) = T(N / size - 1, y);
    }
}


void
border_exchange ( void )
{
    // TODO 7: Communicate border values
    
    // Send my leftmost element to my neighbor on the left
    MPI_Send (
        &(T(1, 1)), 1, MPI_DOUBLE, left_neighbor, 0, MPI_COMM_WORLD
    );
    // Receive my right border element from my right neighbor's leftmost
    MPI_Recv (
        &(T(N / size + 1, M + 1)), 1, MPI_DOUBLE,
        right_neighbor, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE
    );
    // Send my rightmost element to my neighbor on the right
    MPI_Send (
        &(T(N / size, M)), 1, MPI_DOUBLE, right_neighbor, 0, MPI_COMM_WORLD
    );
    // Receive my left border element from my left neighbor's rightmost
    MPI_Recv (
        &(T(0, 0)), 1, MPI_DOUBLE,
        left_neighbor, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE
    );
}


void
domain_init ( void )
{
    // TODO 3: Allocate space for each process' sub-grids
    // and initialize data for the sub-grids
    real_t
        temperature,
        diffusivity;

    left_neighbor = ( rank + size - 1 ) % size;
    right_neighbor = ( rank + size + 1 ) % size;

    len = (N / size + 2) * (M + 2);
    temp[0] = malloc ( len * sizeof(real_t) );
    temp[1] = malloc ( len * sizeof(real_t) );
    thermal_diffusivity = malloc ( len * sizeof(real_t) );

    dt = 0.1;
    dx = 0.1;

    sub_grid_offset = rank * N / size;

    for ( int_t x = sub_grid_offset + 1; x <= N / size + sub_grid_offset; x++ )
    {
        for ( int_t y = 1; y <= M; y++ )
        {
            temperature = 30 + 30 * sin((x + y) / 20.0);
            diffusivity = 0.05 + (30 + 30 * sin((N - x + y) / 20.0)) / 605.0;

            int x_local = x - sub_grid_offset;

            T(x_local,y) = temperature;
            T_next(x_local,y) = temperature;

            THERMAL_DIFFUSIVITY(x_local,y) = diffusivity;
        }
    }

}


void
domain_save ( int_t iteration )
{
    int_t index = iteration / snapshot_frequency;
    char filename[256];
    memset ( filename, 0, 256*sizeof(char) );
    sprintf ( filename, "data/%.5ld.bin", (long) index );

    MPI_File out;
    MPI_File_open(MPI_COMM_WORLD,
                  filename,
                  MPI_MODE_CREATE | MPI_MODE_WRONLY,
                  MPI_INFO_NULL,
                  &out
    );
    if ( ! out ) {
        fprintf(stderr, "Failed to open file: %s\n", filename);
        exit(1);
    }

    MPI_Offset offset = rank * len * sizeof(real_t);
    MPI_File_write_at_all( out, offset, temp[0],
                         len, MPI_DOUBLE , MPI_STATUS_IGNORE ) ;

    MPI_File_close ( &out );
}


void
domain_finalize ( void )
{
    free ( temp[0] );
    free ( temp[1] );
    free ( thermal_diffusivity );
}
