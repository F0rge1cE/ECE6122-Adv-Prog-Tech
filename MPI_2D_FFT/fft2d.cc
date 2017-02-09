// Distributed two-dimensional Discrete FFT transform
// Xueyang Xu
// ECE8893 Project 1


#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <signal.h>
#include <math.h>
#include <mpi.h>

#include "Complex.h"
#include "InputImage.h"

using namespace std;

void Transform1D(Complex* h, int w, Complex* H);
void Transpose(Complex* h, int width, int height);
void inv_Transform1D(Complex* h, int w, Complex* H);
void inv_Transform2D(const char* inputFN);

void Transform2D(const char* inputFN) { 
  // Do the 2D transform here.
  // 1) Use the InputImage object to read in the Tower.txt file and
  //    find the width/height of the input image.
  // 2) Use MPI to find how many CPUs in total, and which one
  //    this process is
  // 3) Allocate an array of Complex object of sufficient size to
  //    hold the 2d DFT results (size is width * height)
  // 4) Obtain a pointer to the Complex 1d array of input data
  // 5) Do the individual 1D transforms on the rows assigned to your CPU
  // 6) Send the resultant transformed values to the appropriate
  //    other processors for the next phase.
  // 6a) To send and receive columns, you might need a separate
  //     Complex array of the correct size.
  // 7) Receive messages from other processes to collect your columns
  // 8) When all columns received, do the 1D transforms on the columns
  // 9) Send final answers to CPU 0 (unless you are CPU 0)
  //   9a) If you are CPU 0, collect all values from other processors
  //       and print out with SaveImageData().
    
  InputImage image(inputFN);  // Create the helper object for reading the image
  // Step (1) in the comments is the line above.
  // Your code here, steps 2-9
    
    int width = image.GetWidth();
    int height = image.GetHeight();
    int nCPUs, myRank, rc;
    int nRows = height; // height is the same as the number of rows

    MPI_Comm_size(MPI_COMM_WORLD, &nCPUs);  //step (2)
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    printf ("Number of tasks= %d My rank= %d\n", nCPUs, myRank);

    int rows_for_each = nRows / nCPUs;
    int my_Start = rows_for_each * myRank;  //starting row for each CPU
    int current_index;                      //indicate the data index in Input_Data

    Complex* Input_Data = image.GetImageData(); //Raw data
    Complex* Output_1D = new Complex[width * height]; //saving 1D-transformed data
    
   
    for (int i = my_Start; i < rows_for_each + my_Start; ++i){
        Complex* h = new Complex[width];
        current_index = i * width;  //point to the begining of each row

        for (int j = 0; j < width; ++j){
            h[j] = Input_Data[current_index + j]; //extract 1-row-data from input
        }

        Complex* H = new Complex[width];   //transformed data, only 1 row
        Transform1D(h, width, H); //do the transform for one row

        for (int k = 0; k < width ; ++k){
            Output_1D[current_index + k] = H[k]; //save the 1-D transform result
        }

        delete h;
        delete H;

    }

    Complex* Receive_1D = new Complex[width * height]; //Complete result receieved by Rank0

    if (myRank != 0){  //other CPUs send data to Rank0
        rc = MPI_Send(&Output_1D[myRank*width*rows_for_each], width*rows_for_each*sizeof(Complex), 
                      MPI_CHAR, 0, 0, MPI_COMM_WORLD);

          if (rc != MPI_SUCCESS)
            {
              cout << "Rank " << myRank
                   << " send failed, rc " << rc << endl;
              MPI_Finalize();
              
            }
          
    }
    else{
      for (int p=0; p<width*rows_for_each; ++p){
          Receive_1D[p] = Output_1D[p];         //Rank0 save its own results
      }
      
    }

    

    if (myRank == 0){    //Rank0 receives from other CPUs

      for (int other_rank = 1; other_rank < nCPUs; ++other_rank){
            MPI_Status status;
            rc = MPI_Recv(&Receive_1D[width*other_rank*rows_for_each], width*rows_for_each*sizeof(Complex), 
                          MPI_CHAR, other_rank, 0, MPI_COMM_WORLD, &status);

            if (rc != MPI_SUCCESS)
              {
                  cout << "Rank " << myRank
                      << " recv failed, rc " << rc << endl;
                  MPI_Finalize();
                  
              }
      } 
     
    }

    Complex* Before_2D = new Complex[width*rows_for_each]; //Used to save data from Rank0

                    //save result after 1-D transform, use Rank0
    if(myRank == 0){         
      image.SaveImageData("MyAfter1D.txt", Receive_1D, width, height);

      // two ways to do the 2D-Transform:
      // (1) Transpose first, repeat the 1D-Trans, then transpose again
      // (2) send by columns(every "width" element), collect and save by columns
      // way (2) needs a for-loop to send and receive, not efficient?
         
      Transpose(Receive_1D, width, height); //Do the transpose, save to Receive_1D

      //Rank0 can save transposed data locally
      for (int t = 0; t< width * rows_for_each; ++t){ 
        Before_2D[t] = Receive_1D[t];
      }

      for (int dest_CPU=1; dest_CPU<nCPUs;++dest_CPU){
        //Send columns from Rank0 to other Ranks

        rc = MPI_Send(&Receive_1D[dest_CPU*width*rows_for_each], width*rows_for_each*sizeof(Complex), 
                      MPI_CHAR, dest_CPU, 0, MPI_COMM_WORLD);

        if (rc != MPI_SUCCESS)
        {
            cout << "Rank " << myRank
                 << " send failed, rc " << rc << endl;
            MPI_Finalize();
           
        }
      }
    }   

    if (myRank !=0){
      MPI_Status status;
      rc = MPI_Recv(&Before_2D[0], width*rows_for_each*sizeof(Complex), 
                          MPI_CHAR, 0, 0, MPI_COMM_WORLD, &status);
      if (rc != MPI_SUCCESS)
        {
           cout << "Rank " << myRank
                << " recv failed, rc " << rc << endl;
           MPI_Finalize();
           
        }
    }

    Complex* After_2D = new Complex[width*rows_for_each];
       
      for(int t1 = 0 ; t1< rows_for_each ; ++t1){ //not robust, width may not equal to height
        Complex* Input_column = new Complex[width];
        int current_index_col = width * t1; //current start point

        for(int t2= 0; t2 < width; ++t2){
          Input_column[t2] = Before_2D[current_index_col+t2];
        }
        Complex* res_one_col_2D = new Complex[width];   //Used to save the result of transform
        
        Transform1D(Input_column, width, res_one_col_2D);

        for (int t3 = 0; t3 < width; ++t3)
        { //save the temp result to After_2D
          After_2D[current_index_col+t3] = res_one_col_2D[t3];
        }
        delete Input_column;
        delete res_one_col_2D;
      }

    // Rank0 collects data from all other Ranks
      Complex* Receive_2D = new Complex[width*height];
      if(myRank != 0){
        rc = MPI_Send(&After_2D[0], width*rows_for_each*sizeof(Complex), 
                      MPI_CHAR, 0, 0, MPI_COMM_WORLD);

        if (rc != MPI_SUCCESS)
          {
            cout << "Rank " << myRank
                   << " send failed, rc " << rc << endl;
            MPI_Finalize();
            
          }
      }

      else{
        //Complex* Receive_2D = new Complex[width*height];
        MPI_Status status;

        for (int i1 = 0; i1 < width*rows_for_each; ++i1){   
            //Save Rank0's own data
          Receive_2D[i1] = After_2D[i1];
        }


        for (int rec_ind = 1; rec_ind < nCPUs; ++rec_ind){
         rc = MPI_Recv(&Receive_2D[width*rec_ind*rows_for_each], width*rows_for_each*sizeof(Complex), 
                           MPI_CHAR, rec_ind, 0, MPI_COMM_WORLD, &status);

         if (rc != MPI_SUCCESS)
           {
              cout << "Rank " << myRank
                   << " recv failed, rc " << rc << endl;
              MPI_Finalize();
              
           }
        }  
      }

      if (myRank == 0){
        //Do the final transpose to get the result
        Transpose(Receive_2D, width, height);

        image.SaveImageData("MyAfter2D.txt",Receive_2D,width,height);

      }

    
    /////////////////////////////////////////////////////////
    // Redo all steps above inversely:
    // 1. transpose first
    // 2. send, process and collect
    // 3. transpose
    // 4. send, process and collect
    // 5. save
      Complex* send_2D_inv = new Complex[width*rows_for_each]; //used as recieve buffers

      if (myRank == 0){
        Transpose(Receive_2D, width, height);

        for (int index_inv = 0; index_inv < rows_for_each*width; ++index_inv){
          send_2D_inv[index_inv] = Receive_2D[index_inv];
        }

        for (int cpu_inv = 1; cpu_inv < nCPUs; ++cpu_inv)
        {
          rc = MPI_Send(&Receive_2D[cpu_inv*width*rows_for_each], width*rows_for_each*sizeof(Complex), 
                      MPI_CHAR, cpu_inv, 0, MPI_COMM_WORLD);

          if (rc != MPI_SUCCESS)
          {
            cout << "Rank " << myRank
                 << " send failed, rc " << rc << endl;
            MPI_Finalize();
           
          }
        }

      }
      else{
        MPI_Status status;
        rc = MPI_Recv(&send_2D_inv[0],rows_for_each*width*sizeof(Complex),MPI_CHAR,0,0,MPI_COMM_WORLD, &status);
        if (rc != MPI_SUCCESS)
        {
          cout << "Rank " << myRank
              << " recv failed, rc " << rc << endl;
          MPI_Finalize();
        }
      }

      
      Complex* result_2D_inv = new Complex[width * rows_for_each];  //used as results of inverse transform

      for(int x = 0; x < rows_for_each; ++x){  //extract one row each time loop
        int current_index_inv = width * x;
        Complex* H_inv = new Complex[width];

        for (int y = 0; y < width ;++y){
          H_inv[y] = send_2D_inv[y+current_index_inv];
        }
        Complex* h_inv = new Complex[width];
        inv_Transform1D(H_inv,width,h_inv);

        for(int z = 0; z < width; ++z){
          result_2D_inv[current_index_inv+z] = h_inv[z];
        }
        delete H_inv;
        delete h_inv;
      }


// //////////////////
// for (int i = my_Start; i < rows_for_each + my_Start; ++i){
//         Complex* h = new Complex[width];
//         current_index = i * width;  //point to the begining of each row

//         for (int j = 0; j < width; ++j){
//             h[j] = Input_Data[current_index + j]; //extract 1-row-data from input
//         }

//         Complex* H = new Complex[width];   //transformed data, only 1 row
//         Transform1D(h, width, H); //do the transform for one row

//         for (int k = 0; k < width ; ++k){
//             Output_1D[current_index + k] = H[k]; //save the 1-D transform result
//         }

//         delete h;
//         delete H;

//     }
// //////////////////













      Complex* Receive_2D_inv = new Complex[width*height]; //data received by thr Rank0

      if(myRank != 0){
        rc = MPI_Send(&result_2D_inv[0], width*rows_for_each*sizeof(Complex), 
                      MPI_CHAR, 0, 0, MPI_COMM_WORLD);
        if (rc != MPI_SUCCESS)
          {
            cout << "Rank " << myRank
                  << " send failed, rc " << rc << endl;
            MPI_Finalize();
          }
      }
      else{
        for (int self_save_inv = 0; self_save_inv < width*rows_for_each; ++self_save_inv){
          Receive_2D_inv[self_save_inv] = result_2D_inv[self_save_inv];
        }

        MPI_Status status;
        for (int rec_indi_inv = 1; rec_indi_inv < nCPUs; ++rec_indi_inv){
          rc = MPI_Recv(&Receive_2D_inv[width*rec_indi_inv*rows_for_each], width*rows_for_each*sizeof(Complex), 
                           MPI_CHAR, rec_indi_inv, 0, MPI_COMM_WORLD, &status);

         if (rc != MPI_SUCCESS)
           {
              cout << "Rank " << myRank
                   << " recv failed, rc " << rc << endl;
              MPI_Finalize();
              
           }
         }
      }
    
    if(myRank == 0){
      Transpose(Receive_2D_inv,width,height);
    }
    //repeat actions like above, process the rows
    Complex* send_1D_inv = new Complex[width*rows_for_each];

    if (myRank == 0){

      for(int x1=0; x1<width*rows_for_each;++x1){
        send_1D_inv[x1] = Receive_2D_inv[x1];
      }
      for(int cpu_inv1=1; cpu_inv1<nCPUs; ++cpu_inv1){
        rc = MPI_Send(&Receive_2D_inv[cpu_inv1*width*rows_for_each], width*rows_for_each*sizeof(Complex), 
                      MPI_CHAR, cpu_inv1, 0, MPI_COMM_WORLD);

          if (rc != MPI_SUCCESS)
          {
            cout << "Rank " << myRank
                 << " send failed, rc " << rc << endl;
            MPI_Finalize();
           
          }
      }
    }
    else{
      MPI_Status status;
        rc = MPI_Recv(&send_1D_inv[0],rows_for_each*width*sizeof(Complex),MPI_CHAR,0,0,MPI_COMM_WORLD, &status);
        if (rc != MPI_SUCCESS)
        {
          cout << "Rank " << myRank
              << " recv failed, rc " << rc << endl;
          MPI_Finalize();
        }
    }
      
    Complex* result_1D_inv = new Complex[width * rows_for_each];

    for(int x1 = 0; x1 < rows_for_each; ++x1){  //extract one row each time loop
      
      int current_index_inv1 = width * x1;

      Complex* H_inv1 = new Complex[width];

      for (int y1 = 0; y1 < width ;++y1){
        H_inv1[y1] = send_1D_inv[y1+current_index_inv1];
      }
      Complex* h_inv1 = new Complex[width];

      inv_Transform1D(H_inv1,width,h_inv1);

      for(int z1 = 0; z1 < width; ++z1){
        result_1D_inv[current_index_inv1+z1] = h_inv1[z1];
      }
      delete h_inv1;
      delete H_inv1;
    }    

    Complex* Receive_1D_inv = new Complex[width * height];

    if(myRank != 0){
        rc = MPI_Send(&result_1D_inv[0], width*rows_for_each*sizeof(Complex), 
                      MPI_CHAR, 0, 0, MPI_COMM_WORLD);
        if (rc != MPI_SUCCESS)
          {
            cout << "Rank " << myRank
                  << " send failed, rc " << rc << endl;
            MPI_Finalize();
          }
      }
      else{
        for (int self_save_inv1 = 0; self_save_inv1 < width*rows_for_each; ++self_save_inv1){
          Receive_1D_inv[self_save_inv1] = result_1D_inv[self_save_inv1];
        }

        MPI_Status status;
        for (int rec_indi_inv1 = 1; rec_indi_inv1 < nCPUs; ++rec_indi_inv1){
          rc = MPI_Recv(&Receive_1D_inv[width*rec_indi_inv1*rows_for_each], width*rows_for_each*sizeof(Complex), 
                           MPI_CHAR, rec_indi_inv1, 0, MPI_COMM_WORLD, &status);
          if (rc != MPI_SUCCESS)
           {
              cout << "Rank " << myRank
                   << " recv failed, rc " << rc << endl;
              MPI_Finalize();
              
           }           
        }
      image.SaveImageData("MyAfterInverse.txt",Receive_1D_inv,width,height);
      
      }
}


void Transform1D(Complex* h, int w, Complex* H){
  // Implement a simple 1-d DFT using the double summation equation
  // given in the assignment handout.  h is the time-domain input
  // data, w is the width (N), and H is the output array.
  int n = 0, i = 0;

  for (n = 0; n < w; ++n){      //the algorithm of Fourier Transform 1D (not inverse)
      H[n] = 0.0;               //initialize the H[n]
      
    for (i = 0; i < w; ++i){
      Complex W( cos(2*M_PI*n*i/w), -sin(2*M_PI*n*i/w) ); //use constructor of Complex
      H[n] = H[n] + h[i] * W ;

    } 
  }  
}


void inv_Transform1D(Complex* H, int w, Complex* h){  
  int n = 0, i = 0;

  for (n = 0; n < w; ++n){      //the algorithm of inverseFourier Transform 
      h[n] = 0.0;               
      
    for (i = 0; i < w; ++i){
      Complex W( cos(2*M_PI*n*i/w)/w, sin(2*M_PI*n*i/w)/w );
      h[n] = h[n] + H[i] * W ;

    } 
  }  
}

void Transpose(Complex* h, int width, int height){
  // Implement a matrix transpose function
  // Only applicable when the matrix is a SQUARE MATRIX!!!
  Complex temp(0,0);

  for (int i = 0; i < height; ++i ){
    for (int j = (i + 1); j < width; ++j){
        temp = h[i*width+j];
        h[i*width+j] = h[j*width+i];
        h[j*width+i] = temp;
    }
  }
}


// void inv_Transform2D(const char* inputFN){
//   InputImage image(inputFN);  // Create the helper object for reading the image
//   // Step (1) in the comments is the line above.
//   // Your code here, steps 2-9
    
//     int width = image.GetWidth();
//     int height = image.GetHeight();
//     int nCPUs, myRank, rc;
//     int nRows = height; // height is the same as the number of rows

//     MPI_Comm_size(MPI_COMM_WORLD, &nCPUs);  //step (2)
//     MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
//     printf ("Number of tasks= %d My rank= %d\n", nCPUs, myRank);

//     int rows_for_each = nRows / nCPUs;
//     int my_Start = rows_for_each * myRank;  //starting row for each CPU
//     int current_index;                      //indicate the data index in Input_Data

//     Complex* Input_Data = image.GetImageData(); //Raw data
//     Complex* Output_1D = new Complex[width * height]; //saving 1D-transformed data

//     //if (myRank == 0){
         
//          //Transpose(Input_Data, width, height);
//     //}
   
//     for (int i = my_Start; i < rows_for_each + my_Start; ++i){
//         Complex* h = new Complex[width];
//         current_index = i * width;  //point to the begining of each row

//         for (int j = 0; j < width; ++j){
//             h[j] = Input_Data[current_index + j]; //extract 1-row-data from input
//         }

//         Complex* H = new Complex[width];   //transformed data, only 1 row
//         inv_Transform1D(h, width, H); //do the transform for one row

//         for (int k = 0; k < width ; ++k){
//             Output_1D[current_index + k] = H[k]; //save the 1-D transform result
//         }

//         delete h;
//         delete H;

//     }

//     Complex* Receive_1D = new Complex[width * height]; //Complete result receieved by Rank0

//     if (myRank != 0){  //other CPUs send data to Rank0
//         rc = MPI_Send(&Output_1D[myRank*width*rows_for_each], width*rows_for_each*sizeof(Complex), 
//                       MPI_CHAR, 0, 0, MPI_COMM_WORLD);

//           if (rc != MPI_SUCCESS)
//             {
//               cout << "Rank " << myRank
//                    << " send failed, rc " << rc << endl;
//               MPI_Finalize();
              
//             }
          
//     }
//     else{
//       for (int p=0; p<width*rows_for_each; ++p){
//           Receive_1D[p] = Output_1D[p];         //Rank0 save its own results
//       }
      
//     }

    

//     if (myRank == 0){    //Rank0 receives from other CPUs

//       for (int other_rank = 1; other_rank < nCPUs; ++other_rank){
//             MPI_Status status;
//             rc = MPI_Recv(&Receive_1D[width*other_rank*rows_for_each], width*rows_for_each*sizeof(Complex), 
//                           MPI_CHAR, other_rank, 0, MPI_COMM_WORLD, &status);

//             if (rc != MPI_SUCCESS)
//               {
//                   cout << "Rank " << myRank
//                       << " recv failed, rc " << rc << endl;
//                   MPI_Finalize();
                  
//               }
//       } 
     
//     }

//     Complex* Before_2D = new Complex[width*rows_for_each]; //Used to save data from Rank0

//                     //save result after 1-D transform, use Rank0
//      if(myRank == 0){         
//        image.SaveImageDataReal("MyAfter1DInverse.txt", Receive_1D, width, height);

//       // two ways to do the 2D-Transform:
//       // (1) Transpose first, repeat the 1D-Trans, then transpose again
//       // (2) send by columns(every "width" element), collect and save by columns
//       // way (2) needs a for-loop to send and receive, not efficient?
         
//       Transpose(Receive_1D, width, height); //Do the transpose, save to Receive_1D

//       //Rank0 can save transposed data locally
//       for (int t = 0; t< width * rows_for_each; ++t){ 
//         Before_2D[t] = Receive_1D[t];
//       }

//       for (int dest_CPU=1; dest_CPU<nCPUs;++dest_CPU){
//         //Send columns from Rank0 to other Ranks

//         rc = MPI_Send(&Receive_1D[dest_CPU*width*rows_for_each], width*rows_for_each*sizeof(Complex), 
//                       MPI_CHAR, dest_CPU, 0, MPI_COMM_WORLD);

//         if (rc != MPI_SUCCESS)
//         {
//             cout << "Rank " << myRank
//                  << " send failed, rc " << rc << endl;
//             MPI_Finalize();
           
//         }
//       }
//     }
       

//     if (myRank !=0){
//       MPI_Status status;
//       rc = MPI_Recv(&Before_2D[0], width*rows_for_each*sizeof(Complex), 
//                           MPI_CHAR, 0, 0, MPI_COMM_WORLD, &status);
//       if (rc != MPI_SUCCESS)
//         {
//            cout << "Rank " << myRank
//                 << " recv failed, rc " << rc << endl;
//            MPI_Finalize();
           
//         }
//     }

//     Complex* After_2D = new Complex[width*rows_for_each];
       
//       for(int t1 = 0 ; t1< rows_for_each ; ++t1){ //not robust, width may not equal to height
//         Complex* Input_column = new Complex[width];
//         int current_index_col = width * t1; //current start point

//         for(int t2= 0; t2 < width; ++t2){
//           Input_column[t2] = Before_2D[current_index_col+t2];
//         }
//         Complex* res_one_col_2D = new Complex[width];   //Used to save the result of transform
        
//         inv_Transform1D(Input_column, width, res_one_col_2D);

//         for (int t3 = 0; t3 < width; ++t3)
//         { //save the temp result to After_2D
//           After_2D[current_index_col+t3] = res_one_col_2D[t3];
//         }
//         delete Input_column;
//         delete res_one_col_2D;
//       }

//     // Rank0 collects data from all other Ranks
//       Complex* Receive_2D = new Complex[width*height];
//       if(myRank != 0){
//         rc = MPI_Send(&After_2D[0], width*rows_for_each*sizeof(Complex), 
//                       MPI_CHAR, 0, 0, MPI_COMM_WORLD);

//         if (rc != MPI_SUCCESS)
//           {
//             cout << "Rank " << myRank
//                    << " send failed, rc " << rc << endl;
//             MPI_Finalize();
            
//           }
//       }

//       else{
        
//         MPI_Status status;

//         for (int i1 = 0; i1 < width*rows_for_each; ++i1){   
//             //Save Rank0's own data
//           Receive_2D[i1] = After_2D[i1];
//         }


//         for (int rec_ind = 1; rec_ind < nCPUs; ++rec_ind){
//          rc = MPI_Recv(&Receive_2D[width*rec_ind*rows_for_each], width*rows_for_each*sizeof(Complex), 
//                            MPI_CHAR, rec_ind, 0, MPI_COMM_WORLD, &status);

//          if (rc != MPI_SUCCESS)
//            {
//               cout << "Rank " << myRank
//                    << " recv failed, rc " << rc << endl;
//               MPI_Finalize();
              
//            }
//         }  
//       }

//       if (myRank == 0){
//         //Do the final transpose to get the result
//         Transpose(Receive_2D, width, height);

//         image.SaveImageDataReal("MyAfterInverse.txt",Receive_2D,width,height);

//       }

// }





int main(int argc, char** argv)
{
  int rc;

  string fn("Tower.txt"); // default file name
  if (argc > 1) fn = string(argv[1]);  // if name specified on cmd line
  
  //string fn1("after2d.txt"); // default file name
  //if (argc > 1) fn1 = string(argv[1]);  // if name specified on cmd line


  rc = MPI_Init(&argc, &argv);   // MPI initialization here
  if (rc != MPI_SUCCESS) {
    printf ("Error starting MPI program. Terminating.\n");
    MPI_Abort(MPI_COMM_WORLD, rc);
  }

  Transform2D(fn.c_str()); // Perform the transform.
 
  //inv_Transform2D(fn1.c_str());

  MPI_Finalize(); // Finalize MPI here

}  
  

  
