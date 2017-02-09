// Threaded two-dimensional Discrete FFT transform
// Xueyang Xu
// ECE8893 Project 2


#include <iostream>
#include <string>
#include <math.h>

#include "Complex.h"
#include "InputImage.h"

#include "pthread.h"  // Optional, default set to be included

// You will likely need global variables indicating how
// many threads there are, and a Complex* that points to the
// 2d image being transformed.

// Global variable visible to all threads
Complex* ImageData;
int      ImageWidth;
int      ImageHeight;
int      nThreads = 16; // Total number of threads
int      ActiveCount;   // Counter of active threads
int      N;

// Mutex variables 
pthread_mutex_t exitMutex; // Used by Transform2D for creating and waiting
pthread_cond_t  exitCond;
pthread_mutex_t ActiveCountMutex;


// Barrier variables and array
int Barrier_count;
pthread_mutex_t Barrier_countMutex;
bool* local_sense = new bool[nThreads];
bool global_sense;

void Transpose(Complex* h, int width, int height);
int FetchAndDecrementCount();
using namespace std;

// Function to reverse bits in an unsigned integer
// This assumes there is a global variable N that is the
// number of points in the 1D transform.
unsigned ReverseBits(unsigned v)
{ //  Provided to students
  unsigned n = N; // Size of array (which is even 2 power k value)
  unsigned r = 0; // Return value
   
  for (--n; n > 0; n >>= 1)
    {
      r <<= 1;        // Shift return value
      r |= (v & 0x1); // Merge in next bit
      v >>= 1;        // Shift reversal value
    }
  return r;
}

// GRAD Students implement the following 2 functions.
// Undergrads can use the built-in barriers in pthreads.

// Call MyBarrier_Init once in main
void MyBarrier_Init()// you will likely need some parameters)
{
  pthread_mutex_init(&Barrier_countMutex,0);
  Barrier_count = nThreads;
  for (int i=0; i<nThreads; i++){
    local_sense[i] = true;
  }
  global_sense = true;
}

// Each thread calls MyBarrier after completing the row-wise DFT
void MyBarrier(unsigned long myID) // Again likely need parameters
{
  local_sense[myID] = !local_sense[myID];
  if (FetchAndDecrementCount() == 1){
    Barrier_count = nThreads;           // Reset the count of threads
    global_sense = !global_sense;
  }
  else{
    while(global_sense != local_sense[myID]) // Spin
      {}
  }
}
  
int FetchAndDecrementCount()
{     // Return the value before decrement!
  pthread_mutex_lock(&Barrier_countMutex);
  int mycount = Barrier_count;
  Barrier_count = Barrier_count - 1;
  pthread_mutex_unlock(&Barrier_countMutex);
  return mycount;
}


void Transform1D(Complex* h, int N)
{
  // Implement the efficient Danielson-Lanczos DFT here.
  // "h" is an input/output parameter
  // "N" is the size of the array (assume even power of 2)


  // 1. Reorder the h, in the reverse code;
  Complex* temp = new Complex[N];
  for (int i = 0; i < N; ++i){
    temp[i] = h[i];
  }
  for (int i = 0; i < N ; ++i){
    h[i] = temp[ReverseBits(i)];
  }
  delete temp; // Avoid the memory leaking

  // 2. Precompute Wn[];
  Complex* W = new Complex[N];

  for (int i = 0; i < N/2; ++i){
    W[i] = Complex(cos(2*M_PI*i/N),-sin(2*M_PI*i/N));
    W[i+N/2] = Complex(-1) * W[i];
  }

  // 3. Do the transformation;
  int Current_length = 2; // sub array begins from the interval of 2
  int Num_of_subArray = N/Current_length; // Number of sub_array in h
  while (Current_length <= N){
    Complex* temp = new Complex[N]();
    for( int i = 0; i < Num_of_subArray; ++i ){
      for( int j = 0; j < Current_length; ++j ){
        temp[i*Current_length + j] = temp[i*Current_length + j] + h[i*Current_length + j%(Current_length/2)] + h[i*Current_length + j%(Current_length/2) + Current_length/2] * W[N*j/Current_length]; 
      }
    }
    for( int j = 0; j < N; ++j ){
      h[j] = temp[j];
    }
    Current_length = Current_length * 2;
    Num_of_subArray = Num_of_subArray / 2;
  delete[] temp;
  }
  delete[] W;
  }

void inv_Transform1D(Complex* h, int N)
{
  // Implement the efficient Danielson-Lanczos DFT here.
  // "h" is an input/output parameter
  // "N" is the size of the array (assume even power of 2)


  // 1. Reorder the h, in the reverse code;
  Complex* temp = new Complex[N];
  for (int i = 0; i < N; ++i){
    temp[i] = h[i];
  }
  for (int i = 0; i < N ; ++i){
    h[i] = temp[ReverseBits(i)];
  }
  delete temp; // Avoid the memory leaking

  // 2. Precompute Wn[];
  Complex* W = new Complex[N];

  for (int i = 0; i < N/2; ++i){
    W[i] = Complex(cos(2*M_PI*i/N),sin(2*M_PI*i/N)); // Difference btw transform and inverse
    W[i+N/2] = Complex(-1) * W[i];
  }

  // 3. Do the transformation;
  int Current_length = 2; // sub array begins from the interval of 2
  int Num_of_subArray = N/Current_length; // Number of sub_array in h
  while (Current_length <= N){
  Complex* temp = new Complex[N]();
    for( int i = 0; i < Num_of_subArray; ++i ){
      for( int j = 0; j < Current_length; ++j ){
        temp[i*Current_length + j] = temp[i*Current_length + j] + h[i*Current_length + j%(Current_length/2)] + h[i*Current_length + j%(Current_length/2) + Current_length/2] * W[N*j/Current_length]; 
      }
    }
    for( int j = 0; j < N; ++j ){
      h[j] = temp[j];
    }
    Current_length *= 2;
    Num_of_subArray /= 2;
    delete[] temp;
  }

  for (int k = 0; k < N; ++k){
    h[k] = h[k]/Complex(N);
    //h[k] = h[k] * Complex(1 / double(N));
    // the "/" operator is not defined in the Complex.h
    // Use * 1/N instead
    // Attention! N is in "int" type, 1/N = 0. Need type conversion.

    if( h[k].Mag().real < 1e-7 ) {
      h[k] = Complex(0);
    } // SaveDataReal() is not defined, use this to judge whether it is real
  
  }

  delete[] W;
  }

void* Transform2DTHread(void* v)
{ // This is the thread startign point.  "v" is the thread number
  // Calculate 1d DFT for assigned rows
  // wait for all to complete
  // Calculate 1d DFT for assigned columns
  // Decrement active count and signal main if all complete
  unsigned long myID = (unsigned long)v;
  int Rows_per_thread = ImageHeight / nThreads; 
  // How many rows for each thread. For this time:1024 / 16 = 64 

  // Complex* h = new Complex[ImageWidth];
  // int start_point = myID * Rows_per_thread;
  // for (i = start_point; i < start_point + Rows_per_thread * ImageWidth ; ++i){
  //   h[i] = ImageData[i];
  // } 
  //It seems no need to create a new array to store the raw data.

  //int start_point = myID * Rows_per_thread; 
  //Start point of data in Imagedata for different threads.

  for (int Row_pointer = 0; Row_pointer < Rows_per_thread; ++Row_pointer){
  // Row_pointer indicates which line is being processed
    Transform1D(ImageData + (myID * Rows_per_thread + Row_pointer) * ImageWidth, ImageWidth);
  }
  //MyBarrier(myID);
  // Finish the transform once (for rows or columns)

  pthread_mutex_lock(&ActiveCountMutex); //Judge whether to send the signal
  ActiveCount--;

  printf("%d\n",ActiveCount); //Test, how many active threads now?

  if (ActiveCount == 0){
    pthread_mutex_unlock(&ActiveCountMutex);
    pthread_mutex_lock(&exitMutex);  //Use this lock to confirm the main thread is ready to receive signal, or it will block
    pthread_cond_signal(&exitCond);
    pthread_mutex_unlock(&exitMutex);
  }
  else{
    pthread_mutex_unlock(&ActiveCountMutex);
  }
  MyBarrier(myID);
  return 0;
}


void* inv_Transform2DTHread(void* v)
{ 
  unsigned long myID = (unsigned long)v;
  int Rows_per_thread = ImageHeight / nThreads; 
  // How many rows for each thread. For this time:1024 / 16 = 64 

  for (int Row_pointer = 0; Row_pointer < Rows_per_thread; ++Row_pointer){
  // Row_pointer indicates which line is being processed
    inv_Transform1D(ImageData + (myID * Rows_per_thread + Row_pointer) * ImageWidth, ImageWidth);
  }
  // Finish the transform once (for rows or columns)

  pthread_mutex_lock(&ActiveCountMutex); //Judge whether to send the signal
  ActiveCount--;

  printf("Inv %d\n",ActiveCount); //Test, how many active threads now?

  if (ActiveCount == 0){
    pthread_mutex_unlock(&ActiveCountMutex);
    pthread_mutex_lock(&exitMutex);  //Use this lock to confirm the main thread is ready to receive signal, or it will block
    pthread_cond_signal(&exitCond);
    pthread_mutex_unlock(&exitMutex);
  }
  else{
    pthread_mutex_unlock(&ActiveCountMutex);
  }
  MyBarrier(myID);
  return 0;
}

void Transform2D(const char* inputFN) 
{ // Do the 2D transform here.
  // Create the helper object for reading the image
  InputImage image(inputFN);  
  // Create the global pointer to the image array data
  ImageData = image.GetImageData();
  ImageWidth = image.GetWidth();
  ImageHeight = image.GetHeight();
  N = ImageWidth;

  // Initialize mutex and condition variables 
  pthread_mutex_init(&exitMutex,0);
  pthread_cond_init(&exitCond, 0);
  MyBarrier_Init();// Initialize the barrier

  //////////////////////////////////////////////
  // Start 16 threads
  pthread_mutex_lock(&exitMutex);
  ActiveCount = nThreads; // Initialize the ActiveCount
  for (int i = 0; i < nThreads; ++i)
    {
      pthread_t pt; // pThread variable (output param from create)
      // Third param is the thread starting function
      // Fourth param is passed to the thread starting function
      pthread_create(&pt, 0, Transform2DTHread, (void*)i);
    }
  pthread_cond_wait(&exitCond, &exitMutex); // Wait for all threads complete

  printf("back successfully!\n"); //Test

  pthread_mutex_unlock(&exitMutex);

  image.SaveImageData("MyAfter1D.txt",ImageData,ImageWidth,ImageHeight);


  //Begin to do the transform of columns
  Transpose(ImageData, ImageWidth, ImageHeight);
  MyBarrier_Init();// Initialize the barrier
  pthread_mutex_lock(&exitMutex);
  ActiveCount = nThreads; // Initialize the ActiveCount
  for (int i = 0; i < nThreads; ++i)
    {
      pthread_t pt; // pThread variable (output param from create)
      // Third param is the thread starting function
      // Fourth param is passed to the thread starting function
      pthread_create(&pt, 0, Transform2DTHread, (void*)i);
    }
  pthread_cond_wait(&exitCond, &exitMutex); // Wait for all threads complete
  pthread_mutex_unlock(&exitMutex);

  Transpose(ImageData, ImageWidth, ImageHeight);

  // Write the transformed data 
  image.SaveImageData("MyAfter2D.txt",ImageData,ImageWidth,ImageHeight);


  // ////////////////////////////////////////////
  // // Do the inverse Transform
  // ////////////////////////////////////////////

  // string fn1("MyAfter2D.txt");
  // //if (argc > 1) fn1 = string(argv[1]);

  // InputImage image1(fn1.c_str());  
  // // Create the global pointer to the image array data
  // ImageData = image1.GetImageData();

  pthread_mutex_lock(&exitMutex);
  ActiveCount = nThreads; // Initialize the ActiveCount
  for (int i = 0; i < nThreads; ++i)
    {
      pthread_t pt; // pThread variable (output param from create)
      // Third param is the thread starting function
      // Fourth param is passed to the thread starting function
      pthread_create(&pt, 0, inv_Transform2DTHread, (void*)i);
    }
  pthread_cond_wait(&exitCond, &exitMutex); // Wait for all threads complete

  printf("Inverse back successfully!\n"); //Test

  pthread_mutex_unlock(&exitMutex);

  image.SaveImageData("MyAfter1DInverse.txt",ImageData,ImageWidth,ImageHeight);


  //Begin to do the transform of columns
  Transpose(ImageData, ImageWidth, ImageHeight);
  MyBarrier_Init();// Initialize the barrier
  pthread_mutex_lock(&exitMutex);
  ActiveCount = nThreads; // Initialize the ActiveCount
  for (int i = 0; i < nThreads; ++i)
    {
      pthread_t pt; // pThread variable (output param from create)
      // Third param is the thread starting function
      // Fourth param is passed to the thread starting function
      pthread_create(&pt, 0, inv_Transform2DTHread, (void*)i);
    }
  pthread_cond_wait(&exitCond, &exitMutex); // Wait for all threads complete
  pthread_mutex_unlock(&exitMutex);

  Transpose(ImageData, ImageWidth, ImageHeight);

  // Write the transformed data 
  image.SaveImageDataReal("MyAfterInverse.txt",ImageData,ImageWidth,ImageHeight);

}


void Transpose(Complex* h, int width, int height)
{
  // Implement a matrix transpose function
  // Only applicable when the matrix is a SQUARE MATRIX!!!
  // Same as Transpose function in homework1.
  Complex temp(0,0);

  for (int i = 0; i < height; ++i ){
    for (int j = (i + 1); j < width; ++j){
        temp = h[i*width+j];
        h[i*width+j] = h[j*width+i];
        h[j*width+i] = temp;
    }
  }
}

int main(int argc, char** argv)
{
  string fn("Tower.txt"); // default file name
  if (argc > 1) fn = string(argv[1]);  // if name specified on cmd line
  
  //local_sense = new bool[nThreads];
  

  Transform2D(fn.c_str()); // Perform the transform.
}  
  

  
