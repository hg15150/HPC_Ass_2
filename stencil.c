
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "mpi.h"

#define NROWS 16
#define MASTER 0

// Define output file name
#define OUTPUT_FILE "stencil.pgm"

void stencil(const int nx, const int ny, float * restrict image, float * restrict tmp_image);
void init_image(const int nx, const int ny, float * image, float *  tmp_image);
void output_image(const char * file_name, const int nx, const int ny, float *image);
double wtime(void);
int calc_nrows_from_rank(int rank, int size, int ny);
void stencilMiddle(const int nx, const int ny, float * restrict image, float * restrict tmp_image);
void sendBottom(const int nx, const int ny, float * restrict tmp_image, int rank);
void sendTop(const int nx, const int ny, float * restrict tmp_image, int rank);
void sendMiddle(const int nx, const int ny, float * restrict tmp_image, int rank);

int tag = 0;
MPI_Status status;
int local_nrows = 0;
int local_ncols = 0;
int remainder = 0;

int main(int argc, char *argv[]) {

  // Check usage
  if (argc != 4) {
    fprintf(stderr, "Usage: %s nx ny niters\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  // Initiliase problem dimensions from command line arguments
  int nx = atoi(argv[1]);     //Dimension x
  int ny = atoi(argv[2]);     //Dimension y
  int niters = atoi(argv[3]); //Number of iterations

  //----------------------------------------------------------------------------
  //MPI section here
  //----------------------------------------------------------------------------

  int rank;
  int above;
  int below;
  int size;
  int flag;
  enum bool {FALSE,TRUE};

  //InitMPI
  MPI_Init(&argc, &argv);
  MPI_Initialized(&flag);
  if ( flag != TRUE ) {
    MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
  }
  MPI_Comm_size( MPI_COMM_WORLD, &size );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  //Calculate size of local rows
  local_nrows = calc_nrows_from_rank(rank, size - 1, ny);
  local_ncols = nx;
  remainder = ny % size;

  above = rank-1;
  below = rank+1;

  //Master branch
  if(rank==MASTER){
    printf("Rank: %d, rows: %d, remainder: %d\n", rank, local_nrows, remainder);


    // Allocate the image
    float *image = malloc(sizeof(float)*nx*ny);
    float *tmp_image = malloc(sizeof(float)*nx*ny);

    init_image(nx, ny, image, tmp_image);

    //Distribute
    int process_size = local_ncols*(local_nrows+1);
    int start_loc = 0;

    //First process
    for (int j = 0; j < process_size; j++) {
        MPI_Ssend(&image[start_loc+j], 1, MPI_FLOAT, 1, tag, MPI_COMM_WORLD);
    }
    // printf("Rank %d 1\n", rank);

    //Middle sections
    process_size = local_ncols*(local_nrows+2);
    for (int i = 2; i < (size-1); i++) {
      start_loc = (i-1)*local_nrows*(local_ncols-1);
      for (int j = 0; j < process_size; j++) {
          MPI_Ssend(&image[start_loc+j],1, MPI_FLOAT, i, tag, MPI_COMM_WORLD);
      }
    }
    // printf("Rank %d 2\n", rank);


    //End section
    start_loc = ((size-1)-1)*local_ncols*local_nrows;
    process_size = local_ncols*(local_nrows+1+remainder);

    for (int j = 0; j < process_size; j++) {
        MPI_Ssend(&image[start_loc+j],1, MPI_FLOAT, size-1, tag, MPI_COMM_WORLD);
    }

    // printf("Rank %d 3\n", rank);

    // Call the stencil kernel
    double tic = wtime();

    //Regather
    int image_portion = local_ncols*local_nrows;
    for (int n_rank = 1; n_rank < size-1; n_rank++) {
      for (int k = 0; k < image_portion; k++) {
          MPI_Recv(&image[(n_rank-1)*image_portion + k], 1, MPI_FLOAT, n_rank, tag, MPI_COMM_WORLD, &status);
      }
    }

    // printf("Rank %d 4\n", rank);


    //Last
    image_portion = local_ncols*(local_nrows+remainder);
    for (int k = 0; k < image_portion; k++) {
        MPI_Recv(&image[((size-1)-1)*image_portion + k], 1, MPI_FLOAT, size-1, tag, MPI_COMM_WORLD, &status);
    }


    double toc = wtime();

    // Output
    printf("------------------------------------\n");
    printf(" runtime: %lf s\n", toc-tic);
    printf("------------------------------------\n");

    output_image(OUTPUT_FILE, nx, ny, image);
    free(image);
    printf("Rank %d\n", rank);
  }

  //First process
  else if(rank == 1){
    // Allocate the image
    float *image = malloc(sizeof(float)*nx*(local_nrows+1));
    float *tmp_image = malloc(sizeof(float)*nx*(local_nrows+1));
    printf("Rank %d: local_nrows = %d\n", rank, local_nrows);


    //Get image from master
    for (int i = 0; i < (local_ncols*(local_nrows+1)); i++) {
      MPI_Recv(&image[i], 1, MPI_FLOAT, MASTER, tag, MPI_COMM_WORLD, &status);
    }

    //Iterate through stencil
    for (int i = 0; i < 1; i++) {
      // stencilMiddle(local_nrows, local_ncols, image, tmp_local_image);
      sendTop(local_ncols, local_nrows, tmp_image, rank);
      // stencilMiddle(local_nrows, local_ncols, tmp_local_image, image);
      sendTop(local_ncols, local_nrows, image, rank);
    }

    //Send image portion back to master
    for (int i = 0; i < (local_ncols*(local_nrows)); i++) {
      MPI_Ssend(&image[i],1, MPI_FLOAT, MASTER, tag, MPI_COMM_WORLD);
    }

    printf("Rank %d\n", rank);
  }

  //Last process
  else if(rank == (size - 1)){
    printf("Rank %d: local_nrows = %d\n", rank, local_nrows);
    // Allocate the image
    float *image = malloc(sizeof(float)*local_ncols*(local_nrows+1));
    float *tmp_image = malloc(sizeof(float)*nx*(local_nrows+1));

    for (int i = 0; i < (local_ncols*(local_nrows+1)); i++) {
      MPI_Recv(&image[i], 1, MPI_FLOAT, MASTER, tag, MPI_COMM_WORLD, &status);
    }

    //Iterate through stencil
    for (int i = 0; i < 1; i++) {
      // stencilMiddle(local_nrows, local_ncols, image, tmp_local_image);
      sendBottom(local_ncols, local_nrows, tmp_image, rank);
      // stencilMiddle(local_nrows, local_ncols, tmp_local_image, image);
      sendBottom(local_ncols, local_nrows, image, rank);
    }

    for (int i = local_ncols; i < (local_ncols*(local_nrows+1)); i++) {
      MPI_Ssend(&image[i],1, MPI_FLOAT, MASTER, tag, MPI_COMM_WORLD);
    }

    printf("Rank %d\n", rank);

  }

  //Middle
  else{
    //Local process image
    float *image = malloc(sizeof(float)*local_ncols*(local_nrows+2));
    float *tmp_local_image = malloc(sizeof(float)*local_ncols*(local_nrows+2));

    printf("Rank %d: local_nrows = %d\n", rank, local_nrows);


    for (int i = 0; i < (local_ncols*(local_nrows+2)); i++) {
      MPI_Recv(&image[i], 1, MPI_FLOAT, MASTER, tag, MPI_COMM_WORLD, &status);
    }
    // printf("Rank %d\n", rank);

    //Iterate through stencil
    for (int i = 0; i < 1; i++) {
      // stencilMiddle(local_nrows, local_ncols, image, tmp_local_image);
      sendMiddle(local_ncols, local_nrows, tmp_local_image, rank);
      // stencilMiddle(local_nrows, local_ncols, tmp_local_image, image);
      sendMiddle(local_ncols, local_nrows, image, rank);
    }

    //Send back to master
    for (int i = local_ncols; i < (local_ncols*(local_nrows+1)); i++) {
      MPI_Ssend(&image[i],1, MPI_FLOAT, MASTER, tag, MPI_COMM_WORLD);
    }
    printf("Rank %d\n", rank);

  }

  MPI_Finalize();
}

//START HERE
//const image

// Stencil given image
void stencil(const int nx, const int ny, float * restrict image, float * restrict tmp_image) {
  //Top left
  tmp_image[0] = image[0]*0.6f + (image[nx] + image[1])*0.1f;
  //Top middle
  for (int i = 1; i < (nx-1); i++) {
    tmp_image[i] = image[i]*0.6f + (image[i-1] + image[i+nx] + image[i+1])*0.1f;
  }
  //Top right
  tmp_image[nx-1] = image[nx-1]*0.6f + (image[nx-2] + image[(nx-1)+nx])*0.1f;

  //Middle section HERE
  for (int r = 1; r < ny-1; r++) {
    //Left
    tmp_image[nx*r] = image[nx*r]*0.6f + (image[nx*(r-1)] + image[nx*(r+1)] + image[nx*r+1])*0.1f;

    //Middle
    for (int c = 1; c < nx-1; c++) {
      tmp_image[c+r*ny] = image[c+r*ny]*0.6f + (image[(c-1)+r*nx] + image[(c+1)+r*nx] + image[c+(r-1)*nx] + image[c+(r+1)*nx])*0.1f;
    }

    //Right
    tmp_image[nx*(r+1)-1] = image[nx*(r+1)-1]*0.6f + (image[nx*r-1] + image[nx*(r+1)-2] + image[nx*(r+2)-1])*0.1f;
  }

  //Bottom left
  tmp_image[(ny-1)*nx] = image[nx*(ny-1)]*0.6f + (image[nx*(ny-2)] + image[(nx*(ny-1))+1])*0.1f;

  //Bottom middle
  for (int i = 1; i < (nx-1); i++) {
    tmp_image[(ny-1)*nx+i] = image[(ny-1)*nx+i]*0.6f + (image[nx*(ny-1)+(i-1)] + image[nx*(ny-1)+(i+1)] + image[(nx)*(ny-2)+i] ) * 0.1f;
  }

  //Bottom right
  tmp_image[ny*nx-1] = image[ny*nx-1]*0.6f + (image[ny*nx-2] + image[ny*(nx-1)-1])*0.1f;

  // for (int j = 0; j < ny; ++j) {
  //   for (int i = 0; i < nx; ++i) {
  //     tmp_image[j+i*ny] = image[j+i*ny] * 0.6f;
  //     if (i > 0)    tmp_image[j+i*ny] += image[j  +(i-1)*ny] * 0.1f;
  //     if (i < nx-1) tmp_image[j+i*ny] += image[j  +(i+1)*ny] * 0.1f;
  //     if (j > 0)    tmp_image[j+i*ny] += image[j-1+i*ny] * 0.1f;
  //     if (j < ny-1) tmp_image[j+i*ny] += image[j+1+i*ny] * 0.1f;
  //   }
  // }
}

// Stencil middle process
void stencilMiddle(const int nx, const int ny, float * restrict image, float * restrict tmp_image){
  for (int r = 1; r < ny-1; r++) {
    //Left
    tmp_image[nx*r] = image[nx*r]*0.6f + (image[nx*(r-1)] + image[nx*(r+1)] + image[nx*r+1])*0.1f;

    //Middle
    for (int c = 1; c < nx-1; c++) {
      tmp_image[c+r*ny] = image[c+r*ny]*0.6f + (image[(c-1)+r*nx] + image[(c+1)+r*nx] + image[c+(r-1)*nx] + image[c+(r+1)*nx])*0.1f;
    }

    //Right
    tmp_image[nx*(r+1)-1] = image[nx*(r+1)-1]*0.6f + (image[nx*r-1] + image[nx*(r+1)-2] + image[nx*(r+2)-1])*0.1f;
  }
}

// Top communication BASIC
void sendTop(const int local_ncols, const int local_nrows, float * restrict tmp_image, int rank){
  if(rank%2==0){
    MPI_Ssend(&tmp_image[local_ncols*(local_nrows+1)], local_ncols, MPI_FLOAT, rank+1, tag, MPI_COMM_WORLD); //Send to rank below
    MPI_Recv(&tmp_image[local_ncols*local_nrows], local_ncols, MPI_FLOAT, rank+1, tag, MPI_COMM_WORLD, &status); //Receive from below
  }
  else {
    MPI_Recv(&tmp_image[local_ncols*local_nrows], local_ncols, MPI_FLOAT, rank+1, tag, MPI_COMM_WORLD, &status); //Receive from below
    MPI_Ssend(&tmp_image[local_ncols*(local_nrows-1)], local_ncols, MPI_FLOAT, rank+1, tag, MPI_COMM_WORLD); //Send to rank below
  }
}

// Bottom communication BASIC
void sendBottom(const int local_ncols, const int local_nrows, float * restrict tmp_image, int rank){
  if(rank%2==0){
    //Top lines
    MPI_Ssend(&tmp_image[local_ncols], local_ncols, MPI_FLOAT, rank-1, tag, MPI_COMM_WORLD); //Send top line to rank above
    MPI_Recv(&tmp_image[0], local_ncols, MPI_FLOAT, rank-1, tag, MPI_COMM_WORLD, &status); //Receive from rank above
  }
  else{
    //Top lines
    MPI_Recv(&tmp_image[0], local_ncols, MPI_FLOAT, rank-1, tag, MPI_COMM_WORLD, &status); //Receive from rank above
    MPI_Ssend(&tmp_image[local_ncols], local_ncols, MPI_FLOAT, rank-1, tag, MPI_COMM_WORLD); //Send top line to rank above
  }
}

// Middle communication BASIC
void sendMiddle(const int local_ncols, const int local_nrows, float * restrict tmp_image, int rank){
  if(rank%2==0){ //Even send up first then receive from above
    //Top lines
    MPI_Ssend(&tmp_image[local_ncols], local_ncols, MPI_FLOAT, rank-1, tag, MPI_COMM_WORLD); //Send top line to rank above
    MPI_Recv(&tmp_image[0], local_ncols, MPI_FLOAT, rank-1, tag, MPI_COMM_WORLD, &status); //Receive from rank above

    //Bottom lines
    MPI_Ssend(&tmp_image[local_ncols*local_nrows], local_ncols, MPI_FLOAT, rank+1, tag, MPI_COMM_WORLD); //Send to rank below
    MPI_Recv(&tmp_image[local_ncols*(local_nrows+1)], local_ncols, MPI_FLOAT, rank+1, tag, MPI_COMM_WORLD, &status); //Receive from below
  }
  else{ //Odd receive first from below then send below
    //Top lines
    MPI_Recv(&tmp_image[0], local_ncols, MPI_FLOAT, rank-1, tag, MPI_COMM_WORLD, &status); //Receive from rank above
    MPI_Ssend(&tmp_image[local_ncols], local_ncols, MPI_FLOAT, rank-1, tag, MPI_COMM_WORLD); //Send top line to rank above

    //Bottom lines
    MPI_Recv(&tmp_image[local_ncols*(local_nrows+1)], local_ncols, MPI_FLOAT, rank+1, tag, MPI_COMM_WORLD, &status); //Receive from below
    MPI_Ssend(&tmp_image[local_ncols*local_nrows], local_ncols, MPI_FLOAT, rank+1, tag, MPI_COMM_WORLD); //Send to rank below
  }
}

// Create the input image
void init_image(const int nx, const int ny, float *  image, float *  tmp_image) {
  // Zero everything
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      image[j+i*ny] = 0.0f;
      tmp_image[j+i*ny] = 0.0f;
    }
  }

  // Checkerboard
  for (int j = 0; j < 8; ++j) {
    for (int i = 0; i < 8; ++i) {
      for (int jj = j*ny/8; jj < (j+1)*ny/8; ++jj) {
        for (int ii = i*nx/8; ii < (i+1)*nx/8; ++ii) {
          if ((i+j)%2)
          image[jj+ii*ny] = 100.0f;
        }
      }
    }
  }
}

// Routine to output the image in Netpbm grayscale binary image format
void output_image(const char * file_name, const int nx, const int ny, float *image) {

  // Open output file
  FILE *fp = fopen(file_name, "w");
  if (!fp) {
    fprintf(stderr, "Error: Could not open %s\n", OUTPUT_FILE);
    exit(EXIT_FAILURE);
  }

  // Ouptut image header
  fprintf(fp, "P5 %d %d 255\n", nx, ny);

  // Calculate maximum value of image
  // This is used to rescale the values
  // to a range of 0-255 for output
  float maximum = 0.0f;
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      if (image[j+i*ny] > maximum)
        maximum = image[j+i*ny];
    }
  }

  // Output image, converting to numbers 0-255
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      fputc((char)(255.0*image[j+i*ny]/maximum), fp);
    }
  }

  // Close the file
  fclose(fp);

}

// Get the current time in seconds since the Epoch
double wtime(void) {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + tv.tv_usec*1e-6;
}

// Calculate the number of rows for each process
int calc_nrows_from_rank(int rank, int size, int ny){
  int nrows;

  nrows = ny / size;       /* integer division */
  if ((ny % size) != 0) {  /* if there is a remainder */
    if (rank == size)
      nrows += ny % size;  /* add remainder to last rank */
  }

  return nrows;
}
