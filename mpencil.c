#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

//TODO
//Test changes to input/output
//Experiment with scaling
//OpenMP?

// Define output file name
#define OUTPUT_FILE "stencil.pgm"
// Iterations defined by niters
//NROWS AND NCOLS ideally not hardcoded
//but otherwise has to be in main??
#define NROWS 1024
#define NCOLS 1024
#define NITERS 200
#define MASTER 0

void init_image(const int nx, const int ny, double *  image);
void output_image(const char * file_name, const int nx, const int ny, double *image);
void stencil(int rank, int size, int lRows, int lCols, double **preImage, double **curImage);
double wtime(void);
int calcNcols(int rank, int size);
int leftCol(int rank, int size);

int main(int argc, char *argv[]) {

  int start;
  int end;
  int size;   //num processes
  int rank;   //rank of process
  int tag = 0; //extra info in send recv
  MPI_Status status; //status struct
  int lRows;
  int lCols; //local number of rows and columns
  int rCols; //remote number of columns
  double **preImage; //local image grid at step iter -1
  double **curImage; //local image grid at step iter
  double *sendbuf; //Send buffer
  double *recvbuf; //Receive buffer
  double *printbuf; //Print buffer

  /*
  double buf[NROWS];
  MPI_File printFile;
  MPI_Request request;
  */

  //Init MPI
  //get size and rank
  MPI_Init( &argc, &argv );
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  //Code from example 11 to do left right calculation
  //Fancy if..else code
  int left = (rank == MASTER) ? (rank + size - 1) : (rank - 1);
  int right = (rank + 1) % size;

  //determine local grid
  lRows = NROWS;
  lCols = calcNcols(rank,size);

  // Allocate the image
  //scale poorly since every process will do this
  double *image = malloc(sizeof(double)*NROWS*NCOLS);

  // Set the input image
  init_image(NROWS, NCOLS, image);
  printf("Rank %d has just run init_image\n",rank);

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
  printf("Rank %d has just allocated memory for arrays\n", rank);

  //initialise local grid

  //gonna be errors due to boundary conditions
  for(int i=0;i<lRows;i++){
    for(int j=0;j<lCols + 1;j++){
      //first row is 0s
      if(i == 0){
        curImage[i][j] = 0.0;
      }
      //last row is 0s
      else if(i == lRows-1){
        curImage[i][j] = 0.0;
      }
      //left column of rank 0 is 0s
      else if((rank == 0) && j == 1){
        curImage[i][j] = 0.0;
      }
      //right column of last rank is 0s
      else if((rank == size-1) && j == lCols){
        curImage[i][j] = 0.0;
      }
      else{
        //equivalent to Boundary Mean case
        curImage[i][j] = image[i+(j+leftCol(rank,size))*NROWS];
      }
    }
  }
  printf("Rank %d has just initialised the local grid with image pixels\n",rank);

  //begin timing
  double tic = wtime();

  for(int iter = 0;iter<NITERS;iter++){
    //send left, receive from right
    for(int i = 0; i < lRows; i++){
      sendbuf[i] = curImage[i][1];
    }
    MPI_Sendrecv( sendbuf, lRows, MPI_DOUBLE, left, tag, recvbuf, lRows, MPI_DOUBLE, right, tag, MPI_COMM_WORLD, &status);
    for(int i = 0; i < lRows; i++){
      curImage[i][lCols+1] = recvbuf[i];
    }

    //send right, receive from left
    for(int i = 0; i < lRows; i++){
      sendbuf[i] = curImage[i][lCols];
    }
    MPI_Sendrecv(sendbuf, lRows, MPI_DOUBLE, right, tag, recvbuf, lRows, MPI_DOUBLE, left, tag, MPI_COMM_WORLD, &status);
    for(int i = 0; i < lRows; i++){
      curImage[i][0] = recvbuf[i];
    }

    //copy old solution into preImage grid
    for(int i = 0;i < lRows; i++){
      for(int j = 0; j < lCols + 2; j++){
        preImage[i][j] = curImage[i][j];
      }
    }

    //run stencil here
    stencil(rank,size,lRows,lCols,preImage,curImage);
  }
  //end timing
  double toc = wtime();
  printf("Rank %d has just finished stencilling\n", rank);

  if(rank == MASTER){
    // Output
    printf("------------------------------------\n");
    printf(" runtime: %lf s\n", toc-tic);
    printf("------------------------------------\n");
  }
  //RANK 0 PROCESS REACHES HERE
  //Do the printing
  double* printImage;
  if(rank == 0){
    start = 2;
    end = lCols;
    printImage = ((double*)malloc(sizeof(double)*NROWS*NCOLS));
    printf("Rank %d has just assigned memory for printImage\n",rank);
  }
  else if(rank == size -1){
    start = 1;
    end = lCols - 1;
  }
  else{
    start = 1;
    end = lCols;
  }
  for(int i = 1; i < lRows-1; i++){
    if(rank == 0){
      printf("Rank %d just got into assigning their own printImage segment\n",rank);
      //change from j = 2; j < lCols + 1
      for(int j = 1; j < lCols; j++){
        printf("Assigning printImage[%d][%d] from curImage[%d][%d]\n",i-1,j-1,i,j);
        //breaks here
        printImage[(i-1*NROWS)+(j-1)] = curImage[i][j];


        printf("Successful Assign of printImage");
      }
      printf("Rank %d just finished their own printImage segment\n",rank);
      for(int k = 1; k < size; k++){
        rCols = calcNcols(k, size);
        MPI_Recv(printbuf, rCols + 2, MPI_DOUBLE, k, tag, MPI_COMM_WORLD, &status);
        printf("Rank %d has just received a printbuffer\n", rank);

        for(int j = 1;j < rCols + 1; j++){
          printImage[(i-1)*NROWS + (j-1)+leftCol(k,size)] = printbuf[j];
        }
      }
    }
    else{
      MPI_Send(curImage[i], lCols + 2, MPI_DOUBLE, MASTER, tag, MPI_COMM_WORLD);
    }
  }
  if(rank==MASTER){
    output_image(OUTPUT_FILE, NROWS, NCOLS, printImage);
  }
  printf("Rank %d has just finished totally\n", rank);

  MPI_Finalize();

  free(image);
  for(int i = 0;i < lRows; i++){
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

// Create the input image
void init_image(const int nx, const int ny, double *  image) {
  // Zero everything
  const int max = nx*ny;
  for (int z = 0; z < max; z++){
    image[z] = 0.0;
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
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      if (image[j+ i*nx] > maximum)
        maximum = image[j + i*nx];
    }
  }

  // Output image, converting to numbers 0-255
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      fputc((char)(255.0*image[j+i*nx]/maximum), fp);
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

void stencil(int rank, int size, int lRows, int lCols, double **preImage, double **curImage){
  int start;
  int end;
  for(int i = 1; i < lRows-1; i++){
    if(rank == 0){
      start = 2;
      end = lCols;
    }
    else if(rank == size -1){
      start = 1;
      end = lCols - 1;
    }
    else{
      start = 1;
      end = lCols;
    }
    for(int j = start; j < end + 1; j++){
      //stencil here
      curImage[i][j] = preImage[i][j] * 0.6;
      curImage[i][j] += preImage[i+1][j] * 0.1;
      curImage[i][j] += preImage[i-1][j] * 0.1;
      curImage[i][j] += preImage[i][j+1] * 0.1;
      curImage[i][j] += preImage[i][j-1] * 0.1;
    }
  }
}
