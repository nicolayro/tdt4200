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

#define MPI_RANK_ROOT  ( rank == 0 )

typedef int64_t int_t;
typedef double real_t;

int 
    rank,           // mpi rank
    size,           // number of ranks
    local_N,        // N in local subgrid
    local_M,        // M in local subgrid
    local_x_offset, // Subgrid offset x
    local_y_offset; // Subgrid offset y


// Cartesean topology
enum Direction { LEFT, RIGHT, DOWN, UP }; // Direction enum for readability
MPI_Comm comm_cart; // Cartesean communicator
int dims[2];        // Dimensions of system [y, x]
int coords[2];      // Position of self in grid system [y, x]
int neighbours[4];  // List of neighbours

int_t
    M,
    N,
    max_iteration,
    snapshot_frequency;

real_t
    *temp[2] = { NULL, NULL },
    *thermal_diffusivity,
    dt;

#define T(x,y)                      temp[0][(y) * (local_N + 2) + (x)]
#define T_next(x,y)                 temp[1][((y) * (local_N + 2) + (x))]
#define THERMAL_DIFFUSIVITY(x,y)    thermal_diffusivity[(y) * (local_N + 2) + (x)]

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
    MPI_Init(&argc, &argv);

    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &size );

    if (MPI_RANK_ROOT) 
    {
        OPTIONS *options = parse_args( argc, argv );

        if ( !options )
        {
            fprintf( stderr, "Argument parsing failed\n" );
            exit(1);
        }

        M = options->M;
        N = options->N;
        max_iteration = options->max_iteration;
        snapshot_frequency = options->snapshot_frequency;
    }

    // Broadcast command line arguments to other ranks.
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&max_iteration, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&snapshot_frequency, 1, MPI_INT, 0, MPI_COMM_WORLD);

    domain_init();

    struct timeval t_start, t_end;
    gettimeofday ( &t_start, NULL );

    for ( int_t iteration = 0; iteration <= max_iteration; iteration++ )
    {

        border_exchange();

        boundary_condition();

        time_step();

        if ( iteration % snapshot_frequency == 0 )
        {
            printf (
                "Iteration %ld of %ld (%.2lf%% complete)\n",
                iteration,
                max_iteration,
                100.0 * (real_t) iteration / (real_t) max_iteration
            );

            domain_save ( iteration );
        }

        swap( &temp[0], &temp[1] );
    }

    gettimeofday ( &t_end, NULL );
    printf ( "Total elapsed time: %lf seconds\n",
            WALLTIME(t_end) - WALLTIME(t_start)
            );


    domain_finalize();
    MPI_Finalize();

    exit ( EXIT_SUCCESS );
}


void
time_step ( void )
{
    real_t c, t, b, l, r, K, new_value;

    // Only looping over local subgrid
    for ( int_t y = 1; y <= local_M; y++ )
    {
        for ( int_t x = 1; x <= local_N; x++ )
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
    // Boundary conditions for cartesean topology
    
    // Bottom
    if (coords[0] == 0) {
        for ( int_t x = 1; x <= local_N; x++ ) {
            T(x, 0) = T(x, 2);
        }
    }

    // Top
    if (coords[0] == dims[0] - 1) {
        for ( int_t x = 1; x <= local_N; x++ ) {
            T(x, local_M+1) = T(x, local_M-1);
        }
    }

    // Left
    if (coords[1] == 0) {
        for (int_t y = 1; y <= local_M; y++) {
            T(0, y) = T(2, y);
        }
    }

    // Right
    if (coords[1] == dims[1] - 1) {
        for (int_t y = 1; y <= local_M; y++) {
            T(local_N+1, y) = T(local_N-1, y);
        }
    }
}


void 
border_exchange( void ) {
    // Setup column and row datatypes for easier border exchange
    MPI_Datatype column, row;
    MPI_Type_vector(local_M, 1, local_N + 2, MPI_DOUBLE, &column);
    MPI_Type_vector(local_N, 1, 1, MPI_DOUBLE, &row);

    MPI_Type_commit ( &column );
    MPI_Type_commit ( &row );

    // Bottom
    MPI_Sendrecv(
        // Send bottom
        &T(1, local_M),
        1,
        row,
        neighbours[UP],
        0,
        // Receive top
        &T(1, 0),
        1,
        row,
        neighbours[DOWN],
        0,
        // Communicator and status
        comm_cart,
        MPI_STATUS_IGNORE
    );

    // Up
    MPI_Sendrecv(
        // Send up
        &T(1, 1),
        1,
        row,
        neighbours[DOWN],
        0,
        // Receive down
        &T(1, local_M + 1),
        1,
        row,
        neighbours[UP],
        0,
        // Communicator and status
        comm_cart,
        MPI_STATUS_IGNORE
    );

    // Left
    MPI_Sendrecv(
        // Send left
        &T(1, 1),
        1,
        column,
        neighbours[LEFT],
        0,
        // Receive right
        &T(local_N + 1, 1),
        1,
        column,
        neighbours[RIGHT],
        0,
        // Communicator and status
        comm_cart,
        MPI_STATUS_IGNORE
    );

    // Right
    MPI_Sendrecv(
        // Send right
        &T(local_N, 1),
        1,
        column,
        neighbours[RIGHT],
        0,
        // Receive Left
        &T(0, 1),
        1,
        column,
        neighbours[LEFT],
        0,
        // Communicator and status
        comm_cart,
        MPI_STATUS_IGNORE
    );
}

void
domain_init ( void )
{
    int periods[2] = {0, 0};

    // Dimensions
    MPI_Dims_create(size, 2, dims);
    // Create cartesian grid with new communicator
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &comm_cart);

    // Get coordinates 
    MPI_Cart_coords(comm_cart, rank, 2, coords);

    // Local values and offsets
    local_N = N / dims[1]; // x
    local_M = M / dims[0]; // y
    local_x_offset = coords[1] * local_N;
    local_y_offset = coords[0] * local_M;

    // Get neighbours;
    // 0: left, 1: right, 2: down, 3: up
    MPI_Cart_shift(comm_cart, 1, 1, &neighbours[LEFT], &neighbours[RIGHT]);
    MPI_Cart_shift(comm_cart, 0, 1, &neighbours[DOWN], &neighbours[UP]);

    // Allocate memory for calculations - only fit local subgrid + halo
    temp[0] = malloc ( (local_M+2)*(local_N+2) * sizeof(real_t) );
    temp[1] = malloc ( (local_M+2)*(local_N+2) * sizeof(real_t) );
    thermal_diffusivity = malloc ( (local_M+2)*(local_N+2) * sizeof(real_t) );

    dt = 0.1;

    for ( int_t y = 1; y <= local_M; y++ )
    {
        for ( int_t x = 1; x <= local_N; x++ )
        {
            // We offset x and y values for calculations
            real_t x_1 = x + local_x_offset;
            real_t y_1 = y + local_y_offset;

            real_t temperature = 30 + 30 * sin((x_1 + y_1) / 20.0);
            real_t diffusivity = 0.05 + (30 + 30 * sin((N - x_1 + y_1) / 20.0)) / 605.0;

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

    // Data type for contain data. By using subarray we can easily remove
    // the halo of the local subgrid.
    MPI_Datatype data;
    int size2[2] = { local_M + 2, local_N + 2};
    int grid2[2] = { local_M, local_N };
    int start2[2] = { 1, 1 };
    MPI_Type_create_subarray(2, size2, grid2, start2,
                             MPI_ORDER_C, MPI_DOUBLE, &data);
    MPI_Type_commit(&data);

    // Setup datatype to represent subgrid location in global world.
    MPI_Datatype area;
    int size[2] = { M, N };
    int grid[2] = { local_M, local_N };
    int start[2] = { local_y_offset, local_x_offset };
    MPI_Type_create_subarray(2, size, grid, start, 
                             MPI_ORDER_C, MPI_DOUBLE, &area);
    MPI_Type_commit(&area);

    MPI_File out;
    MPI_File_open(
        comm_cart, 
        filename,
        MPI_MODE_CREATE | MPI_MODE_WRONLY,
        MPI_INFO_NULL,
        &out
    );

    if ( ! out ) {
        fprintf(stderr, "Failed to open file: %s\n", filename);
        exit(1);
    }

    // Restrict writing from this rank to the sub-area we defined above
    MPI_File_set_view (
        out, 0, MPI_DOUBLE, area, "native", MPI_INFO_NULL
    );

    MPI_File_write_all(
        out,
        temp[0],
        1,
        data,
        MPI_STATUS_IGNORE
    );

    MPI_File_close(&out);
}


void
domain_finalize ( void )
{
    free ( temp[0] );
    free ( temp[1] );
    free ( thermal_diffusivity );
}
