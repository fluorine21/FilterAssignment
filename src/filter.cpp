#include <iostream>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <cstdlib>
#include "bmplib.h"

using namespace std;

//============================Add function prototypes here======================

void convolve(unsigned char out[][SIZE][RGB], unsigned char in[][SIZE][RGB], int N, double kernel[][11]);
void dummy(unsigned char out[][SIZE][RGB], unsigned char in[][SIZE][RGB]);
void sobel(unsigned char out[][SIZE][RGB], unsigned char in[][SIZE][RGB]);
void gaussian(double kernel[][11], int N, double sigma);
void gaussian_filter(unsigned char out[][SIZE][RGB], unsigned char in[][SIZE][RGB], int N, double sigma);
void unsharp(unsigned char out[][SIZE][RGB], unsigned char in[][SIZE][RGB], int N, double sigma, double alpha);

//============================Do not change code in main()======================

#ifndef AUTOTEST

int main(int argc, char* argv[])
{
   //First check argc
  if(argc < 3)
    {
      //we need at least ./filter <input file> <filter name> to continue
      cout << "usage: ./filter <input file> <filter name> <filter parameters>";
      cout << " <output file name>" << endl;
      return -1;
    }
   //then check to see if we can open the input file
   unsigned char input[SIZE][SIZE][RGB];
   unsigned char output[SIZE][SIZE][RGB];
   char* outfile;
   int N;
   double sigma, alpha;
   double kernel[11][11];

   // read file contents into input array
   int status = readRGBBMP(argv[1], input);
   if(status != 0)
   {
      cout << "unable to open " << argv[1] << " for input." << endl;
      return -1;
   }
   //Input file is good, now look at next argument
   if( strcmp("sobel", argv[2]) == 0)
   {
     sobel(output, input);
     outfile = argv[3];
   }
   else if( strcmp("blur", argv[2]) == 0)
   {
     if(argc < 6)
       {
	 cout << "not enough arguments for blur." << endl;
	 return -1;
       }
     N = atoi(argv[3]);
     sigma = atof(argv[4]);
     outfile = argv[5];
     gaussian_filter(output, input, N, sigma);
   }
   else if( strcmp("unsharp", argv[2]) == 0)
   {
     if(argc < 7)
       {
	 cout << "not enough arguments for unsharp." << endl;
	 return -1;
       }
     N = atoi(argv[3]);
     sigma = atof(argv[4]);
     alpha = atof(argv[5]);
     outfile = argv[6];
     unsharp(output, input, N, sigma, alpha);

   }
   else if( strcmp("dummy", argv[2]) == 0)
   {
     //do dummy
     dummy(output, input);
     outfile = argv[3];
   }
   else
   {
      cout << "unknown filter type." << endl;
      return -1;
   }

   if(writeRGBBMP(outfile, output) != 0)
   {
      cout << "error writing file " << argv[3] << endl;
   }

}   

#endif 

//=========================End Do not change code in main()=====================


// Creates an identity kernel (dummy kernel) that will simply
// copy input to output image and applies it via convolve()
//
// ** This function is complete and need not be changed.
// Use this as an example of how to create a kernel array, fill it in
// appropriately and then use it in a call to convolve.
void dummy(unsigned char out[][SIZE][RGB], unsigned char in[][SIZE][RGB])
{
	double k[11][11];
	for (int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			k[i][j] = 0;
		}
	}
	k[1][1] = 1;
	convolve(out, in, 3, k);
}


// Convolves an input image with an NxN kernel to produce the output kernel
// You will need to complete this function by following the 
//  instructions in the comments
void convolve(unsigned char out[][SIZE][RGB], unsigned char in[][SIZE][RGB], 
		int N, double kernel[][11])
{

	int padded[SIZE+10][SIZE+10][RGB];  // Use for input image with appropriate
	// padding
	int temp[SIZE][SIZE][RGB];          // Use for the unclamped output pixel
	// values then copy from temp to out,
	// applying clamping


	//first set all of padded to 0 (black)
	for (int i = 0; i < (SIZE + 10); i++) {
		for (int j = 0; j < (SIZE + 10); j++) {
			for (int k = 0; k < RGB; k++) {
				padded[i][j][k] = 0;
			}
		}
	}


	//now copy input into padded to appropriate locations
	for (int i = 5; i < SIZE+5; i++) {
		for (int j = 5; j < SIZE+5; j++) {
			for (int k = 0; k < RGB; k++) {
				padded[i][j][k] = in[i-5][j-5][k];
			}
		}
	}


	//initialize temp pixels to 0 (black)
	for (int i = 0; i < SIZE; i++) {
		for (int j = 0; j < SIZE; j++) {
			for (int k = 0; k < RGB; k++) {
				temp[i][j][k] = 0;
			}
		}
	}


	//now perform convolve (using convolution equation on each pixel of the
	// actual image) placing the results in temp (i.e. unclamped result)
	//Here we give you the structure of the convolve for-loops, you need
	//to figure out the loop limits
	for (int y = 5; y < SIZE+5; y++) {
		for (int x = 5; x < SIZE+5; x++) {
			for (int k = 0; k < RGB; k++) {
				for (int i = -(N/2); i <= (N/2) ; i++) {
					for (int j = -(N/2); j <= (N/2); j++) {
						temp[x-5][y-5][k] += (padded[x+i][y+j][k] * kernel[(N/2)+i][(N/2)+j]);
					}
				}
			}
		}
	}


	//now clamp and copy to output
	// You may need to cast to avoid warnings from the compiler:
	for (int i = 0; i < SIZE; i++) {
		for (int j = 0; j < SIZE; j++) {
			for (int k = 0; k < RGB; k++) {
			   if(temp[i][j][k] > 0 && temp[i][j][k] < 265){
					out[i][j][k] = (unsigned char)temp[i][j][k];
				}
				else if(temp[i][j][k] > 255){
				   out[i][j][k] = 255;
				}
				else{
				   out[i][j][k] = 0;
				}
			}
		}
	}
	// (i.e. out[i][j][k] = (unsigned char) temp[i][j][k];)




}

// You will need to complete this function by following the 
//  instructions in the comments
void sobel(unsigned char out[][SIZE][RGB], unsigned char in[][SIZE][RGB])
{
	double k[11][11];
	double s_h1[3][3] = { {-1, 0, 1},
			{-2, 0, 2},
			{-1, 0, 1} };
	double s_h2[3][3] = { {1, 0, -1},
			{2, 0, -2},
			{1, 0, -1} };
	unsigned char h1_soble[SIZE][SIZE][RGB]; //hold intemediate images
	unsigned char h2_soble[SIZE][SIZE][RGB];

	for (int i = 0; i < 11; i++)
	{
		for(int j = 0; j < 11; j++)
		{
			k[i][j] = 0;
		}
	}

	// Copy in 1st 3x3 horizontal sobel kernel (i.e. copy s_h1 into k)
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			k[i][j] = s_h1[i][j];
		}
	}


	// Call convolve to apply horizontal sobel kernel with result in h1_soble
	convolve(h1_soble, in, 3, k);


	// Copy in 2nd 3x3 horizontal sobel kernel (i.e. copy s_h2 into k)
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			k[i][j] = s_h2[i][j];
		}
	}


	// Call convolve to apply horizontal sobel kernel with result in h2_soble
	convolve(h2_soble, in, 3, k);




	// Add the two results (applying clamping) to produce the final output

	int temp[SIZE][SIZE][RGB];

	for (int i = 0; i < SIZE; i++) {
		for (int j = 0; j < SIZE; j++) {
			for (int l = 0; l < RGB; l++) {
				temp[i][j][l] = (h1_soble[i][j][l] + h2_soble[i][j][l]);
			}
		}
	}

	for (int i = 0; i < SIZE; i++) {
		for (int j = 0; j < SIZE; j++) {
			for (int l = 0; l < RGB; l++) {
			   if(temp[i][j][l] > 0 && temp[i][j][l] < 265){
					out[i][j][l] = (unsigned char)temp[i][j][l];
				}
				else if(temp[i][j][l] > 255){
				   out[i][j][l] = 255;
				}
				else{
				   out[i][j][l] = 0;
				}
			}
		}
	}


}

//Generates the Kernel to be used, stores it in the 2D array passed to it (hopefully)
void gaussian(double kernel[][11], int N, double sigma){

	//A = 1 so we're ignoring it
	double r, s = 2.0 * sigma * sigma;

	// sum will be used to normalize kernel
	double sum = 0.0;

	// generate NxN kernel
	for (int x = (int)(-N/2); x <= (int)(N/2); x++)
	{
		for(int y = (int)(-N/2); y <= (int)(N/2); y++)
		{
			//r stores our exponent
			r = (x*x + y*y)/s;
			//Final calculation is made here and Gaussian indexes are converted to C++ indexes
			kernel[x + (int)(N/2)][y + (int)(N/2)] = exp(r*-1.0);
			sum += kernel[x + (int)(N/2)][y + (int)(N/2)];
		}
	}

	// normalize the Kernel
	for(int i = 0; i < N; ++i){
		for(int j = 0; j < N; ++j){
			kernel[i][j] /= sum;
		}
	}

   cout << endl;
	//Prints the final kernel to the screen
	for(int i = 0; i < N; i++){
		for(int j = 0; j < N; j++){
			cout << kernel[i][j] << " ";
		}
		cout << endl;
	}
	sum = 0;
	for(int i = 0; i < N; i++){
		for(int j = 0; j < N; j++){
			sum += kernel[i][j];
		}
	}

	cout << "Gaussian Sum = " << sum << endl;
	cout << endl << "This is the end, my only friend, the end" << endl;
}
void gaussian_filter(unsigned char out[][SIZE][RGB], unsigned char in[][SIZE][RGB], int N, double sigma){
	//Passes a blank kernel to the gaussian method in order to generate a filter kernel
	double kernel[11][11];
	gaussian(kernel, N, sigma);
	convolve(out, in, N, kernel);


}
void unsharp(unsigned char out[][SIZE][RGB], unsigned char in[][SIZE][RGB], int N, double sigma, double alpha){

	//Stores a blurred image in "out"
	gaussian_filter(out, in, N, sigma);

   int D[SIZE][SIZE][RGB];
   int S[SIZE][SIZE][RGB];

	//Subtracts the blurred image from the original input image and stores it in temp
	for(int i = 0; i < SIZE; i++){
		for(int j = 0; j < SIZE; j++){
			for(int k = 0; k < RGB; k++){
				D[i][j][k] = in[i][j][k] - out[i][j][k];
			}
		}
	}

	for(int i = 0; i < SIZE; i++){
		for(int j = 0; j < SIZE; j++){
			for(int k = 0; k < RGB; k++){
				S[i][j][k] = in[i][j][k] + (alpha*(D[i][j][k]));
			}
		}
	}



	//copies values of in into out
	for (int i = 0; i < SIZE; i++) {
		for (int j = 0; j < SIZE; j++) {
			for (int k = 0; k < RGB; k++) {
			   if(S[i][j][k] > 0 && S[i][j][k] < 265){
					out[i][j][k] = (unsigned char)S[i][j][k];
				}
				else if(S[i][j][k] > 255){
				   out[i][j][k] = 255;
				}
				else{
				   out[i][j][k] = 0;
				}
			}
		}
	}

}
