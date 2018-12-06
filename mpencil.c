
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

// Define output file name
#define OUTPUT_FILE "stencil.pgm"

void stencil(const int nx, const int ny, double *  image, double *  tmp_image);
void init_image(const int nx, const int ny, double *  image, double *  tmp_image);
void output_image(const char * file_name, const int nx, const int ny, double *image);
double wtime(void);

int main(int argc, char *argv[]) {

  // Check usage
  if (argc != 4) {
    fprintf(stderr, "Usage: %s nx ny niters\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  // Initiliase problem dimensions from command line arguments
  int nx = atoi(argv[1]);
  int ny = atoi(argv[2]);
  int niters = atoi(argv[3]);

  // Allocate the image
  double *image = malloc(sizeof(double)*nx*ny);
  double *tmp_image = malloc(sizeof(double)*nx*ny);

  // Set the input image
  init_image(nx, ny, image, tmp_image);

  // Call the stencil kernel
  double tic = wtime();
  for (int t = 0; t < niters; ++t) {
    stencil(nx, ny, image, tmp_image);
    stencil(nx, ny, tmp_image, image);
  }
  double toc = wtime();


  // Output
  printf("------------------------------------\n");
  printf(" runtime: %lf s\n", toc-tic);
  printf("------------------------------------\n");

  output_image(OUTPUT_FILE, nx, ny, image);
  free(image);
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
void init_image(const int nx, const int ny, double *  image, double *  tmp_image) {
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
