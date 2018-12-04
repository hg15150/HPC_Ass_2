
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
  int tag = 0;
  int local_nrows;
  int local_ncols;
  int size;

  //InitMPI
  MPI_Init(&argc, &argv);
  MPI_Comm_size( MPI_COMM_WORLD, &size );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  //Calculate size of local rows
  local_ncols = nx;
  local_nrows = calc_nrows_from_rank(rank, size - 1, ny);

  //Master branch
  if(rank==MASTER){

    // Allocate the image
    float *image = malloc(sizeof(float)*nx*ny);
    float *tmp_image = malloc(sizeof(float)*nx*ny);

    init_image(nx, ny, image, tmp_image);

    // Call the stencil kernel
    double tic = wtime();

    //Distribute
    for (int i = 1; i < (rank-1); i++) {
      for (int j = 0; j < (local_ncols*(local_nrows+2))-1; j++) {


      }


    }

    //Regather

    double toc = wtime();

  }

  else if(rank == 1){
    //Set ranks of the processes above and below
    below = rank + 1;

  }

  else if(rank == (size - 1)){


  }

  else{
    //Local process image
    float *local_image = malloc(sizeof(float)*local_ncols*(local_nrows+2));
    float *tmp_local_image = malloc(sizeof(float)*local_ncols*(local_nrows+2));


  }


  for (int t = 0; t < niters; ++t) {
    stencil(nx, ny, image, tmp_image);
    stencil(nx, ny, tmp_image, image);
  }

  MPI_Finalize();



  // Output
  printf("------------------------------------\n");
  printf(" runtime: %lf s\n", toc-tic);
  printf("------------------------------------\n");

  output_image(OUTPUT_FILE, nx, ny, image);
  free(image);
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
    if (rank == size - 1)
      nrows += NROWS % size;  /* add remainder to last rank */
  }

  return nrows;
}
