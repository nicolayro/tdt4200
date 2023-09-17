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
    rank;
int_t *local_sizes = NULL;

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

    N = options->N;
    M = options->M;
    max_iteration = options->max_iteration;
    snapshot_frequency = options->snapshot_frequency;

    // Task 2
    int sendData[M][N], receiveData[M][N];
    if (rank == 0) {
        for (int i = 0; i < size; i++) {
            MPI_Send(sendData, N, MPI_INT, 1, 0, MPI_COMM_WORLD);
            MPI_Recv(receiveData, N, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    } else {
        MPI_Recv(receiveData, N, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
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
                iteration,
                max_iteration,
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

    for ( int_t x = 1; x <= N; x++ )
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
    for ( int_t x = 1; x <= N; x++ )
    {
        T(x, 0) = T(x, 2);
        T(x, M+1) = T(x, M-1);
    }

    for ( int_t y = 1; y <= M; y++ )
    {
        T(0, y) = T(2, y);
        T(N+1, y) = T(N-1, y);
    }
}


void
border_exchange ( void )
{
    // TODO 7: Communicate border values
}


void
domain_init ( void )
{
    // TODO 3: Allocate space for each process' sub-grids
    // and initialize data for the sub-grids
    local_sizes = malloc(N * sizeof(int_t));

    for (int_t i = 0; i < size; i++) {
        local_sizes[i] = (int_t) N / size;
    }

    real_t
        temperature,
        diffusivity;

    temp[0] = malloc ( (N+2)*(M+2) * sizeof(real_t) );
    temp[1] = malloc ( (N+2)*(M+2) * sizeof(real_t) );
    thermal_diffusivity = malloc ( (N+2)*(M+2) * sizeof(real_t) );

    dt = 0.1;
    dx = 0.1;

    for ( int_t x = 1; x <= N; x++ )
    {
        for ( int_t y = 1; y <= M; y++ )
        {
            temperature = 30 + 30 * sin((x + y) / 20.0);
            diffusivity = 0.05 + (30 + 30 * sin((N - x + y) / 20.0)) / 605.0;
            T(x,y) = temperature;
            T_next(x,y) = temperature;

            THERMAL_DIFFUSIVITY(x,y) = diffusivity;
        }
    }
}


void
domain_save ( int_t iteration )
{
    int_t index = iteration / snapshot_frequency;
    char filename[256];
    memset ( filename, 0, 256*sizeof(char) );
    sprintf ( filename, "data/%.5ld.bin", index );

    FILE *out = fopen ( filename, "wb" );
    if ( ! out ) {
        fprintf(stderr, "Failed to open file: %s\n", filename);
        exit(1);
    }
    fwrite( temp[0], sizeof(real_t), (N+2)*(M+2), out );
    fclose ( out );
}


void
domain_finalize ( void )
{
    free ( temp[0] );
    free ( temp[1] );
    free ( thermal_diffusivity );
}
