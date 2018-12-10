#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

// Define output file name
#define OUTPUT_FILE "stencil.pgm"
// Iterations defined by niters
//NROWS AND NCOLS ideally not hardcoded
//but otherwise has to be in main??
#define NROWS 1024
#define NCOLS 1024
#define NITERS 200
#define MASTER 0

void stencil(const int nx, const int ny, double *  image, double *  tmp_image);
void init_image(const int nx, const int ny, double *  image);
void output_image(const char * file_name, const int nx, const int ny, double *image);
double wtime(void);
int calcNcols(int rank, int size);
int leftCol(int rank);

int main(int argc, char *argv[]) {

  int i,j;    //row and column indices
  int k;      //rank iterator
  int start,end; //start and end columns - RANK dependent
  int iter;   //timestep iterator
  int rank;   //rank of process
  int left,right; //left right rank trackers
  int size;   //num processes
  int tag = 0; //extra info in send recv
  MPI_Status status; //status struct
  int lRows,lCols; //local number of rows and columns
  int rCols; //remote number of columns
  double **preImage; //local image grid at step iter -1
  double **curImage; //local image grid at step iter
  double *sendbuf; //Send buffer
  double *recvbuf; //Receive buffer
  double *printbuf; //Print buffer


  // Check usage
  if (argc != 4) {
    fprintf(stderr, "Usage: %s nx ny niters\n", argv[0]);
    exit(EXIT_FAILURE);
  }


  //Init MPI
  //get size and rank
  MPI_Init( &argc, &argv );
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  //Code from example 11 to do left right calculation
  //Fancy if..else code
  left = (rank == MASTER) ? (rank + size - 1) : (rank - 1);
  right = (rank + 1) % size;

  //determine local grid
  lRows = NROWS;
  lCols = calcNcols(rank,size);

  // Allocate the image
  //scale poorly since every process will do this
  double *image = malloc(sizeof(double)*NROWS*NCOLS);

  // Set the input image
  init_image(NROWS, NCOLS, image);

  //allocate space for local grid
  //two columns added for HALO
  preImage = (double**)malloc(sizeof(double*) * lRows);
  for(int i=0;i<lRows;i++){
    preImage[i] = (double*)malloc(sizeof(double) * (lCols+2));
  }
  curImage = (double**)malloc(sizeof(double*) * lRows);
  for(int i=0;i<lRows;i++){
    curImage[i] = (double*)malloc(sizeof(double) * (lCols+2));
  }
  sendbuf = (double*)malloc(sizeof(double) * lRows);
  recvbuf = (double*)malloc(sizeof(double) * lRows);

  //last rank has most columns in this format
  rCols = calcNcols(size-1,size);
  printbuf = (double*)malloc(sizeof(double) * (rCols + 2));

  //initialise local grid

  //gonna be errors due to boundary conditions
  for(int i=0;i<lRows;i++){
    for(int j=0;j<lCols + 1;j++){
      if(i = 0){
        curImage[i][j] = 0.0;
      }
      elseif(i == lRows-1){
        curImage[i][j] = 0.0;
      }
      elseif((rank == 0) && j == 1){
        curImage[i][j] = 0.0;
      }
      elseif((rank == size-1) && j == lCols){
        curImage[i][j] = 0.0;
      }
      else{
        //equivalent to Boundary Mean case
        curImage[i][j] = image[i][j+leftCol(rank)]
      }
    }
  }

  //begin timing
  double tic = wtime();

  for(iter = 0;iter<NITERS;iter++){
    //send left, receive from right
    for(i = 0; i < lRows; i++){
      sendbuf[i] = curImage[i][1];
    }
    MPI_Sendrecv( sendbuf, lRows, MPI_DOUBLE, left, tag, recvbuf, lRows, MPI_DOUBLE, right, tag, MPI_COMM_WORLD, &status);
    for(i = 0; i < lRows; i++){
      curImage[i][lCols+1] = recvbuf[i];
    }

    //send right, receive from left
    for(i = 0; i < lRows; i++){
      sendbuf[i] = curImage[i][lCols];
    }
    MPI_Sendrecv(sendbuf, lRows, MPI_DOUBLE, right, rag, recvbuf, lRows, MPI_DOUBLE, left, tag, MPI_COMM_WORLD, &status);
    for(i = 0; i < lRows; i++){
      curImage[i][0] = recvbuf[i];
    }

    //copy old solution into preImage grid
    for(i = 0;i < lRows; i++){
      for(j = 0; j < lCols + 2; j++){
        preImage[i][j] = curImage[i][j];
      }
    }

    //run stencil here
    for(i = 1; i < lRows-1; i++){
      if(rank == 0){
        start = 2;
        end = lCols;
      }
      elseif(rank == size -1){
        start = 1;
        end = lCols - 1;
      }
      else{
        start = 1;
        end = lCols;
      }
      for(j = start; j < end + 1; j++){
        //stencil here
        curImage[i][j] = preImage[i][j] * 0.6;
        curImage[i][j] += preImage[i+1][j] * 0.1;
        curImage[i][j] += preImage[i-1][j] * 0.1;
        curImage[i][j] += preImage[i][j+1] * 0.1;
        curImage[i][j] += preImage[i][j-1] * 0.1;
      }
    }
  }
  //end timing
  double toc = wtime();

  /*
  for (int t = 0; t < niters; ++t) {
    stencil(nx, ny, image, tmp_image);
    stencil(nx, ny, tmp_image, image);
  }
  */

  // Output
  printf("------------------------------------\n");
  printf(" runtime: %lf s\n", toc-tic);
  printf("------------------------------------\n");

  //Do the printing
  for(i = 1; i < lRows-1; i++){
    if(rank == 0){
      start = 2;
      end = lCols;
    }
    elseif(rank == size -1){
      start = 1;
      end = lCols - 1;
    }
    else{
      start = 1;
      end = lCols;
    }
    if(rank == 0){
      for(j = 2; j < lCols + 1; j++){
        image[i][j] = curImage[i][j]
      }
      for(k = 1; k < size; k++){
        rCols = calcNcols(k, size);
        MPI_Recv(printbuf, rCols + 2, MPI_DOUBLE, k, tag, MPI_COMM_WORLD, &status);
        for(j = 1;j < rCols + 1; j++){
          image[i][(j-1)+leftCol] = printbuf[j];
        }
      }
    }
    else{
      MPI_Send(curImage[i], lCols + 2, MPI_DOUBLE, MASTER, tag, MPI_COMM_WORLD, &status);
    }
  }

  output_image(OUTPUT_FILE, NROWS, NCOLS, image);

  MPI_Finalize();

  free(image);
  for(i = 0;i < lRows; i++){
    free(preImage[i]);
    free(curImage[i]);
  }
  free(preImage);
  free(curImage);
  free(sendbuf);
  free(recvbuf);
  free(printbuf);

  return EXIT_SUCCESS;
}

//called 200 times
void stencil(const int nx, const int ny, double * restrict image, double * restrict tmp_image) {
  const int max = nx*ny;
  //Remove conditionals by precomputing
  //corners
  tmp_image[0] = image[0] * 0.6;
  tmp_image[0] += image[1] * 0.1;
  tmp_image[0] += image[nx] * 0.1;

  tmp_image[nx-1] = image[nx-1] * 0.6;
  tmp_image[nx-1] += image[nx-2] *0.1;
  tmp_image[nx-1] += image[2*nx -1] *0.1;

  tmp_image[max-nx] = image[max-nx] * 0.6;
  tmp_image[max-nx] += image[max-(2*nx)] * 0.1;
  tmp_image[max-nx] += image[max-nx+1] * 0.1;

  tmp_image[max-1] = image[max-1] * 0.6;
  tmp_image[max-1] += image[max-2] * 0.1;
  tmp_image[max-1] += image[max-nx-1] * 0.1;

  //top row
  for (int t = 1; t < (nx-1); t++){
    tmp_image[t] = image[t] * 0.6;
    tmp_image[t] += image[t-1] * 0.1;
    tmp_image[t] += image[t+1] * 0.1;
    tmp_image[t] += image[t+nx] * 0.1;
  }

  //left column
  for (int l = nx; l <= (max-nx); l+=nx){
    tmp_image[l] = image[l] *0.6;
    tmp_image[l] += image[l-nx] *0.1;
    tmp_image[l] += image[l+1] *0.1;
    tmp_image[l] += image[l+nx] *0.1;
  }

  //right column
  for (int r = (nx-1); r < (max-nx); r+=nx){
    tmp_image[r] = image[r] *0.6;
    tmp_image[r] = image[r-nx] *0.1;
    tmp_image[r] = image[r-1] *0.1;
    tmp_image[r] = image[r+nx] *0.1;
  }

  //bottom row
  for (int b = (max-nx)+1; b < (max-1); b++){
    tmp_image[b] = image[b] *0.6;
    tmp_image[b] += image[b-1] *0.1;
    tmp_image[b] += image[b+1] *0.1;
    tmp_image[b] += image[b-nx] *0.1;
  }

  //main body square
  //need to figure out cycle
  /*
  for (int z = nx+1; z<(max-nx)-1; z++){
    tmp_image[z] = image[z] *0.6;
    tmp_image[z] += image[z-nx] *0.1;
    tmp_image[z] += image[z-1] *0.1;
    tmp_image[z] += image[z+1] *0.1;
    tmp_image[z] += image[z+nx] *0.1;
    if(z%nx == (nx-2)){
      z+=2;
    }
  }
  */
  for(int i = 1; i<ny-1; i++){
    for(int j = 1; j<nx-1; j++){
      tmp_image[j+i*nx] = image[j+i*nx] *0.6;
      tmp_image[j+i*nx] += image[j+((i-1)*nx)] *0.1;
      tmp_image[j+i*nx] += image[(j-1)+i*nx] *0.1;
      tmp_image[j+i*nx] += image[(j+1)+i*nx] *0.1;
      tmp_image[j+i*nx] += image[j+((i+1)*nx)] *0.1;
    }
  }

  /*
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      tmp_image[j+i*nx] = image[j+i*nx] * 0.6;
      if (i > 0)    tmp_image[j+i*nx] += image[j  +(i-1)*nx] * 0.1;
      if (i < ny-1) tmp_image[j+i*nx] += image[j  +(i+1)*nx] * 0.1;
      if (j > 0)    tmp_image[j+i*nx] += image[j-1+i*nx] * 0.1;
      if (j < nx-1) tmp_image[j+i*nx] += image[j+1+i*nx] * 0.1;
    }
  }
}*/
/*
void stencil(const int nx, const int ny, double * image, double * tmp_image) {
  const int max = nx*ny;
  const int modnix = (max%nx)-1;
  for(int z = 0; z <= max; z++){
    tmp_image[z] = image[z] * 0.6;
    if (z-nx > 0)    tmp_image[z] += image[z-nx] * 0.1;
    if (z+nx < max) tmp_image[z] += image[z+nx] * 0.1;
    if (z % nx > 0)    tmp_image[z] += image[z-1] * 0.1;
    if (z < modnix) tmp_image[z] += image[z+1] * 0.1;
  }
}
*/
}
// Create the input image
void init_image(const int nx, const int ny, double *  image) {
  // Zero everything
  const int max = nx*ny;
  for (int z = 0; z < max; z++){
    image[z] = 0.0;
    tmp_image[z] = 0.0;
  }

  // Checkerboard
  for (int j = 0; j < 8; ++j) {
    for (int i = 0; i < 8; ++i) {
      for (int jj = j*ny/8; jj < (j+1)*ny/8; ++jj) {
        for (int ii = i*nx/8; ii < (i+1)*nx/8; ++ii) {
          if ((i+j)%2)
          image[jj+ii*ny] = 100.0;
        }
      }
    }
  }
}

// Routine to output the image in Netpbm grayscale binary image format
void output_image(const char * file_name, const int nx, const int ny, double *image) {

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
  double maximum = 0.0;
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

int calcNcols(int rank, int size){
  int ncols;

  ncols = NCOLS / size;       /* integer division */
  if ((NCOLS % size) != 0) {  /* if there is a remainder */
    if (rank == size - 1)
      ncols += NCOLS % size;  /* add remainder to last rank */
  }

  return ncols;
}

int leftCol(int rank, int size){
  int l;
  l = rank * (NCOLS/size);
  return l;
}
