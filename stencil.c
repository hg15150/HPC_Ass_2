
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
void stencilMiddle(float * restrict image, float * tmp_image);
void stencilTop(float * restrict image, float * tmp_image);
void stencilBottom(float * restrict image, float * tmp_image);
void sendBottom(float * restrict tmp_image, int rank);
void sendTop(float * restrict tmp_image, int rank);
void sendMiddle(float * restrict tmp_image, int rank);

int tag = 0;
MPI_Status status;
int R;
int n_rows;
int P;
int Rem;
int rem_rows;

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
  // niters = 100;

  //----------------------------------------------------------------------------
  //MPI section here
  //----------------------------------------------------------------------------

  int rank;
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

  R = nx;
  n_rows = ny/(size);
  P = R*n_rows;
  Rem = (ny%(size))*nx;
  rem_rows = ny%(size);


  //Master branch
  if(rank==MASTER){
    // Allocate the image
    float *image = malloc(sizeof(float)*nx*ny);
    float *tmp_image = malloc(sizeof(float)*nx*ny);

    init_image(nx, ny, image, tmp_image);

    //Distribute
    // int process_size = P+R;

    // //First process
    // for (int j = 0; j < process_size; j++) {
    //     MPI_Ssend(&image[j], 1, MPI_FLOAT, 1, tag, MPI_COMM_WORLD);
    // }

    //Middle sections
    int process_size = P+2*R;
    int start_loc = P-R;
    for (int i = 1; i < (size-1); i++) {
      start_loc = (i*P)-R;
      for (int j = 0; j < process_size; j++) {
          MPI_Ssend(&image[start_loc+j],1, MPI_FLOAT, i, tag, MPI_COMM_WORLD);
      }
    }

    //End section
    start_loc = ((size-1)*P)-R;
    process_size = P+R+Rem;

    for (int j = 0; j < process_size; j++) {
        MPI_Ssend(&image[start_loc+j],1, MPI_FLOAT, size-1, tag, MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double tic = wtime();
    //Iterate through stencil
    for (int i = 0; i < niters; i++) {
      stencilTop(image, tmp_image);
      sendTop(tmp_image, rank);
      stencilTop(tmp_image, image);
      sendTop(image, rank);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double toc = wtime();

    //Regather
    for (int n_rank = 1; n_rank < size-1; n_rank++) {
      int start_loc_recv = n_rank*P;
      for (int k = 0; k < P; k++) {
          MPI_Recv(&image[start_loc_recv + k], 1, MPI_FLOAT, n_rank, tag, MPI_COMM_WORLD, &status);
      }
    }

    //Last
    int image_portion = P+Rem;
    int start_recv = (size-1)*P;
    for (int k = 0; k < image_portion; k++) {
        MPI_Recv(&image[start_recv + k], 1, MPI_FLOAT, size-1, tag, MPI_COMM_WORLD, &status);
    }

    // Output
    printf("------------------------------------\n");
    printf(" runtime: %lf s\n", toc-tic);
    printf("------------------------------------\n");

    output_image(OUTPUT_FILE, nx, ny, image);
    free(image);
  }

  /*
  //First process
  /*
  else if(rank == 1){
    // Allocate the image
    float *image = malloc(sizeof(float)*(P+R));
    float *tmp_image = malloc(sizeof(float)*(P+R));

    //Get image from master
    for (int i = 0; i < P+R; i++) {
      MPI_Recv(&image[i], 1, MPI_FLOAT, MASTER, tag, MPI_COMM_WORLD, &status);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    //Iterate through stencil
    for (int i = 0; i < niters; i++) {
      stencilTop(image, tmp_image);
      sendTop(tmp_image, rank);
      stencilTop(tmp_image, image);
      sendTop(image, rank);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    //Send image portion back to master
    for (int i = 0; i < P; i++) {
      MPI_Ssend(&image[i],1, MPI_FLOAT, MASTER, tag, MPI_COMM_WORLD);
    }
  }
  */

  //Last process
  else if(rank == (size - 1)){
    // Allocate the image
    float *image = malloc(sizeof(float)*(P+R+Rem));
    float *tmp_image = malloc(sizeof(float)*(P+R+Rem));

    for (int i = 0; i < (P+R+Rem); i++) {
      MPI_Recv(&image[i], 1, MPI_FLOAT, MASTER, tag, MPI_COMM_WORLD, &status);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    //Iterate through stencil
    for (int i = 0; i < niters; i++) {
      stencilBottom(image, tmp_image);
      sendBottom(tmp_image, rank);
      stencilBottom(tmp_image, image);
      sendBottom(image, rank);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    for (int i = 0; i < P+Rem; i++) {
      MPI_Ssend(&image[R+i],1, MPI_FLOAT, MASTER, tag, MPI_COMM_WORLD);
    }
  }

  //Middle
  else{
    //Local process image
    float *image = malloc(sizeof(float)*(P+2*R));
    float *tmp_local_image = malloc(sizeof(float)*(P+2*R));

    for (int i = 0; i < (P+2*R); i++) {
      MPI_Recv(&image[i], 1, MPI_FLOAT, MASTER, tag, MPI_COMM_WORLD, &status);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    //Iterate through stencil
    for (int i = 0; i < niters; i++) {
      stencilMiddle(image, tmp_local_image);
      sendMiddle(tmp_local_image, rank);
      stencilMiddle(tmp_local_image, image);
      sendMiddle(image, rank);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    //Send back to master
    for (int i = 0; i < P; i++) {
      MPI_Ssend(&image[R+i],1, MPI_FLOAT, MASTER, tag, MPI_COMM_WORLD);
    }
  }

  MPI_Finalize();
  printf("Rank %d FINISHED\n", rank);
}

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

// Stencil bottom
void stencilBottom(float * restrict image, float * restrict tmp_image){
  //Middle section HERE
  for (int r = 1; r < (n_rows+rem_rows); r++) {
    //Left
    tmp_image[R*r] = image[R*r]*0.6f + (image[R*(r-1)] + image[R*(r+1)] + image[R*r+1])*0.1f;
    //Middle
    for (int c = 1; c < R-1; c++) {
      tmp_image[c+r*R] = image[c+r*R]*0.6f + (image[(c-1)+r*R] + image[(c+1)+r*R] + image[c+(r-1)*R] + image[c+(r+1)*R])*0.1f;
    }
    //Right
    tmp_image[R*(r+1)-1] = image[R*(r+1)-1]*0.6f + (image[R*r-1] + image[R*(r+1)-2] + image[R*(r+2)-1])*0.1f;
  }

  //Bottom left
  tmp_image[(n_rows+rem_rows)*R] = image[R*(n_rows+rem_rows)]*0.6f + (image[R*(n_rows+rem_rows-1)] + image[(R*(n_rows+rem_rows))+1])*0.1f;

  //Bottom middle
  for (int i = 1; i < (R-1); i++) {
    tmp_image[(n_rows+rem_rows)*R+i] = image[(n_rows+rem_rows)*R+i]*0.6f + (image[R*(n_rows+rem_rows)+(i-1)] + image[R*(n_rows+rem_rows)+(i+1)] + image[(R)*(n_rows+rem_rows-1)+i] ) * 0.1f;
  }

  //Bottom right
  tmp_image[(n_rows+rem_rows+1)*R-1] = image[(n_rows+rem_rows+1)*R-1]*0.6f + (image[(n_rows+rem_rows+1)*R-2] + image[(n_rows+rem_rows)*R-1])*0.1f;
}

// Stencil top
void stencilTop(float * restrict image, float * restrict tmp_image){
  //Top left
  tmp_image[0] = image[0]*0.6f + (image[R] + image[1])*0.1f;
  //Top middle
  for (int i = 1; i < (R-1); i++) {
    tmp_image[i] = image[i]*0.6f + (image[i-1] + image[i+R] + image[i+1])*0.1f;
  }
  //Top right
  tmp_image[R-1] = image[R-1]*0.6f + (image[R-2] + image[(R-1)+R])*0.1f;

  //Middle section HERE
  for (int r = 1; r < n_rows; r++) {
    //Left
    tmp_image[R*r] = image[R*r]*0.6f + (image[R*(r-1)] + image[R*(r+1)] + image[R*r+1])*0.1f;
    //Middle
    for (int c = 1; c < R-1; c++) {
      tmp_image[c+r*R] = image[c+r*R]*0.6f + (image[(c-1)+r*R] + image[(c+1)+r*R] + image[c+(r-1)*R] + image[c+(r+1)*R])*0.1f;
    }
    //Right
    tmp_image[R*(r+1)-1] = image[R*(r+1)-1]*0.6f + (image[R*r-1] + image[R*(r+1)-2] + image[R*(r+2)-1])*0.1f;
  }
}

// Stencil middle process
void stencilMiddle(float * restrict image, float * restrict tmp_image){
  for (int j = 1; j <= n_rows; j++) {
    //Left
    tmp_image[R*j] = image[R*j]*0.6f + (image[R*(j-1)] + image[R*(j+1)] + image[R*j+1])*0.1f;
    //Middle
    for (int i = 1; i < R-1; i++) {
      tmp_image[i+j*R] = image[i+j*R]*0.6f + (image[(i-1)+j*R] + image[(i+1)+j*R] + image[i+(j-1)*R] + image[i+(j+1)*R])*0.1f;
    }
    //Right
    tmp_image[R*(j+1)-1] = image[R*(j+1)-1]*0.6f + (image[R*j-1] + image[R*(j+1)-2] + image[R*(j+2)-1])*0.1f;
  }
}

// Top communication BASIC
void sendTop(float * restrict tmp_image, int rank){

  if(rank%2==0){
    MPI_Ssend(&tmp_image[P-R], R, MPI_FLOAT, rank+1, tag, MPI_COMM_WORLD); //Send to rank below
    MPI_Recv(&tmp_image[P], R, MPI_FLOAT, rank+1, tag, MPI_COMM_WORLD, &status); //Receive from below
  }
  else {
    MPI_Recv(&tmp_image[P], R, MPI_FLOAT, rank+1, tag, MPI_COMM_WORLD, &status); //Receive from below
    MPI_Ssend(&tmp_image[P-R], R, MPI_FLOAT, rank+1, tag, MPI_COMM_WORLD); //Send to rank below
  }
}

// Bottom communication BASIC
void sendBottom(float * restrict tmp_image, int rank){
  if(rank%2==0){
    MPI_Ssend(&tmp_image[R], R, MPI_FLOAT, rank-1, tag, MPI_COMM_WORLD); //Send to rank below
    MPI_Recv(&tmp_image[0], R, MPI_FLOAT, rank-1, tag, MPI_COMM_WORLD, &status); //Receive from below
  }
  else {
    MPI_Recv(&tmp_image[0], R, MPI_FLOAT, rank-1, tag, MPI_COMM_WORLD, &status); //Receive from below
    MPI_Ssend(&tmp_image[R], R, MPI_FLOAT, rank-1, tag, MPI_COMM_WORLD); //Send to rank below
  }
}

// Middle communication BASIC
void sendMiddle(float * restrict tmp_image, int rank){
  if(rank%2==0){ //Even send up first then receive from above
    //Top lines
    MPI_Ssend(&tmp_image[R], R, MPI_FLOAT, rank-1, tag, MPI_COMM_WORLD); //Send top line to rank above
    MPI_Recv(&tmp_image[0], R, MPI_FLOAT, rank-1, tag, MPI_COMM_WORLD, &status); //Receive from rank above

    //Bottom lines
    MPI_Ssend(&tmp_image[P], R, MPI_FLOAT, rank+1, tag, MPI_COMM_WORLD); //Send to rank below
    MPI_Recv(&tmp_image[P+R], R, MPI_FLOAT, rank+1, tag, MPI_COMM_WORLD, &status); //Receive from below
  }
  else{ //Odd receive first from below then send below
    //Bottom lines
    MPI_Recv(&tmp_image[P+R], R, MPI_FLOAT, rank+1, tag, MPI_COMM_WORLD, &status); //Receive from below
    MPI_Ssend(&tmp_image[P], R, MPI_FLOAT, rank+1, tag, MPI_COMM_WORLD); //Send to rank below

    //Top lines
    MPI_Recv(&tmp_image[0], R, MPI_FLOAT, rank-1, tag, MPI_COMM_WORLD, &status); //Receive from rank above
    MPI_Ssend(&tmp_image[R], R, MPI_FLOAT, rank-1, tag, MPI_COMM_WORLD); //Send top line to rank above
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
  if ((ny % size) != 0) {  /* if there is a rem */
    if (rank == size)
      nrows += ny % size;  /* add rem to last rank */
  }

  return nrows;
}
