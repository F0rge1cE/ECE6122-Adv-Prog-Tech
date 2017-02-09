// Distributed two-dimensional Discrete FFT transform
// XI WU GTID-902849967
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

void Transform1D(Complex* h, int w, Complex* H, int flag);

void Transform2D(const char* inputFN)
{ // Do the 2D transform here.
  // 1) Use the InputImage object to read in the Tower.txt file and
  //    find the width/height of the input image.
  // 2) Use MPI to find how many CPUs in total, and which one
  //    this process is
  // number of rows per CPU
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
   
  // Create the helper object for reading the image
  // Step (1) in the comments is the line above.
  // Your code here, steps 2-9

    InputImage image(inputFN);
    int width=image.GetWidth();
    int height=image.GetHeight();

    int rc;
    int nCPUs=0;
    int rank=0;
    MPI_Comm_size(MPI_COMM_WORLD,&nCPUs);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    printf("Number of tasks= %d My rank= %d\n",nCPUs,rank);
    int RowCPU=height/nCPUs;

    Complex* Res1D = new Complex[width*height];
    Complex* DataIn;
    DataIn=image.GetImageData(); 
  // define current index and row to indicate the position
    int CurrentIndex;
    int CurrentRow = rank*RowCPU;
  // 1D Fourier transform for rows
    for (int i=CurrentRow;i<CurrentRow+RowCPU;i++)
    {
      CurrentIndex = width*i; 
      Complex* h = new Complex[width];
  // extract the row that needs to perform fft
      for (int j=0;j<width;j++)
      {
        h[j] = DataIn[CurrentIndex+j];
      }
      Complex* H = new Complex[width];
      Transform1D(h,width,H,0);
  // save the outputs into a matrix
      for (int k=0;k<width;k++)
      {
        Res1D[CurrentIndex+k] = H[k];
      }
    }
    delete(DataIn);

    Complex* Recv1D = new Complex[width*height];
  // other cpus send to cpu 0.
    if (rank!=0)
    {
        rc = MPI_Send(&Res1D[rank*width*RowCPU],RowCPU*width*sizeof(Complex),MPI_CHAR,0,0,MPI_COMM_WORLD);
    }
    else
    {
  // save cpu 0's own data.
      for (int r=0;r<RowCPU*width;r++)
      {
        Recv1D[r] = Res1D[r];
      }
  // cpu 0 receives from other cpus.
      MPI_Status status;
      for (int cpu=1;cpu<nCPUs;cpu++)
      {
        rc = MPI_Recv(&Recv1D[cpu*width*RowCPU],RowCPU*width*sizeof(Complex),MPI_CHAR,cpu,0,MPI_COMM_WORLD, &status);
        if (rc != MPI_SUCCESS)
        {
          cout << "Rank " << rank
               << " recv failed, rc " << rc << endl;
          MPI_Finalize();
        }
      }
    }

// Send transposed columns to other cpus for 1D transform
    Complex* Send2D = new Complex[width*RowCPU];
    if (rank==0)
    {
      image.SaveImageData("MyAfter1d.txt",Recv1D,width,height);
// Matrix Transpose
      Complex temp;
      for (int t1=0;t1<height;t1++)
      {
        for (int t2=(t1+1);t2<width;t2++)
        {
          temp = Recv1D[width*t1+t2];
          Recv1D[width*t1+t2] = Recv1D[width*t1+t2+(t2-t1)*(width-1)];
          Recv1D[width*t1+t2+(t2-t1)*(width-1)] = temp;
        }
        for (int ind=0;ind<RowCPU*width;ind++)
        {
          Send2D[ind] = Recv1D[ind];
        }
      }  
      for (int cpu0=1;cpu0<nCPUs;cpu0++)
      {
        rc = MPI_Send(&Recv1D[cpu0*width*RowCPU],RowCPU*width*sizeof(Complex),MPI_CHAR,cpu0,0,MPI_COMM_WORLD);
  // MPI send failure
        if (rc != MPI_SUCCESS)
        {
          cout << "Rank " << rank
               << " send failed, rc " << rc << endl;
          MPI_Finalize();
        }
      }
    }
    else
    {
      MPI_Status status;
      rc = MPI_Recv(&Send2D[0],RowCPU*width*sizeof(Complex),MPI_CHAR,0,0,MPI_COMM_WORLD, &status);
      if (rc != MPI_SUCCESS)
      {
        cout << "Rank " << rank
             << " recv failed, rc " << rc << endl;
        MPI_Finalize();
      }
    }

// define current index to indicate the position
    int CurrentIndex1;
// 1D Fourier transform for transposed columns
    Complex* Res2D = new Complex[width*RowCPU];
    for (int i1=0;i1<RowCPU;i1++)
    {
      CurrentIndex1 = width*i1; 
      Complex* h1 = new Complex[width];
  // extract the transposed columns that needs to perform fft
      for (int j1=0;j1<width;j1++)
      {
        h1[j1] = Send2D[CurrentIndex1+j1];
      }
      Complex* H1 = new Complex[width];
      Transform1D(h1,width,H1,0);
  // save the outputs into a matrix
      for (int k1=0;k1<width;k1++)
      {
        Res2D[CurrentIndex1+k1] = H1[k1];
      }
    }

    Complex* Recv2D = new Complex[width*height];
  // other cpus send to cpu 0.
    if (rank!=0)
    {
        rc = MPI_Send(&Res2D[0],RowCPU*width*sizeof(Complex),MPI_CHAR,0,0,MPI_COMM_WORLD);
  // MPI send failure
        if (rc != MPI_SUCCESS)
        {
          cout << "Rank " << rank
               << " send failed, rc " << rc << endl;
          MPI_Finalize();
        }
    }
    else
    {
  // save cpu 0's own data.
      for (int r1=0;r1<RowCPU*width;r1++)
      {
        Recv2D[r1] = Res2D[r1];
      }
  // cpu 0 receives from other cpus.
      MPI_Status status;
      for (int cpu1=1;cpu1<nCPUs;cpu1++)
      {
        rc = MPI_Recv(&Recv2D[cpu1*width*RowCPU],RowCPU*width*sizeof(Complex),MPI_CHAR,cpu1,0,MPI_COMM_WORLD, &status);
        if (rc != MPI_SUCCESS)
        { 
          cout << "Rank " << rank
               << " recv failed, rc " << rc << endl;
          MPI_Finalize();
        }
      }
    }
    if (rank==0)
    {
  // Matrix Transpose
      Complex temp1;
      for (int t3=0;t3<height;t3++)
      {
        for (int t4=(t3+1);t4<width;t4++)
        {
          temp1 = Recv2D[width*t3+t4];
          Recv2D[width*t3+t4] = Recv2D[width*t3+t4+(t4-t3)*(width-1)];
          Recv2D[width*t3+t4+(t4-t3)*(width-1)] = temp1;
        }
      }
      image.SaveImageData("MyAfter2d.txt", Recv2D, width, height);
    }

  

  // Inverse FFT2D
  // Transposed Matrix in cpu 0
    Complex* Send2DInv = new Complex[width*RowCPU];
    if (rank==0)
    {
      Complex tempinv;
      for (int ti1=0;ti1<height;ti1++)
      {
        for (int ti2=(ti1+1);ti2<width;ti2++)
        {
          tempinv = Recv2D[width*ti1+ti2];
          Recv2D[width*ti1+ti2] = Recv2D[width*ti1+ti2+(ti2-ti1)*(width-1)];
          Recv2D[width*ti1+ti2+(ti2-ti1)*(width-1)] = tempinv;
        }
      }
      for (int indi=0;indi<RowCPU*width;indi++)
      {
        Send2DInv[indi] = Recv2D[indi];
      }
      for (int cpui=1;cpui<nCPUs;cpui++)
      {
        rc = MPI_Send(&Recv2D[cpui*width*RowCPU],RowCPU*width*sizeof(Complex),MPI_CHAR,cpui,0,MPI_COMM_WORLD);
  // MPI send failure
        if (rc != MPI_SUCCESS)
        {
          cout << "Rank " << rank
               << " send failed, rc " << rc << endl;
          MPI_Finalize();
        }
      }
    }
    else
    {
      MPI_Status status;
      rc = MPI_Recv(&Send2DInv[0],RowCPU*width*sizeof(Complex),MPI_CHAR,0,0,MPI_COMM_WORLD, &status);
      if (rc != MPI_SUCCESS)
      {
        cout << "Rank " << rank
             << " recv failed, rc " << rc << endl;
        MPI_Finalize();
      }
    }

    int CurrentIndexInv;
    Complex* Res2DInv = new Complex[width*RowCPU];
    for (int ii=0;ii<RowCPU;ii++)
    {
      CurrentIndexInv = width*ii; 
      Complex* Hi = new Complex[width];
  // extract the transposed columns that needs to perform fft
      for (int ji=0;ji<width;ji++)
      {
        Hi[ji] = Send2DInv[CurrentIndexInv+ji];
      }
      Complex* hi = new Complex[width];
      Transform1D(Hi,width,hi,1);
  // save the outputs into a matrix
      for (int ki=0;ki<width;ki++)
      {
        Res2DInv[CurrentIndexInv+ki] = hi[ki];
      }
    }

    Complex* Recv2DInv = new Complex[width*height];
  // other cpus send to cpu 0.
    if (rank!=0)
    {
      rc = MPI_Send(&Res2DInv[0],RowCPU*width*sizeof(Complex),MPI_CHAR,0,0,MPI_COMM_WORLD);
  // MPI send failure
      if (rc != MPI_SUCCESS)
      {
        cout << "Rank " << rank
             << " send failed, rc " << rc << endl;
        MPI_Finalize();
      }
    }
    else
    {
  // save cpu 0's own data.
      for (int ri=0;ri<RowCPU*width;ri++)
      {
        Recv2DInv[ri] = Res2DInv[ri];
      }
  // cpu 0 receives from other cpus.
      MPI_Status status;
      for (int cpui1=1;cpui1<nCPUs;cpui1++)
      {
        rc = MPI_Recv(&Recv2DInv[cpui1*width*RowCPU],RowCPU*width*sizeof(Complex),MPI_CHAR,cpui1,0,MPI_COMM_WORLD, &status);
        if (rc != MPI_SUCCESS)
        {
          cout << "Rank " << rank
               << " recv failed, rc " << rc << endl;
          MPI_Finalize();
        }
      }
    }
    if (rank==0)
    {
  // Matrix Transpose
      Complex tempinv1;
      for (int ti3=0;ti3<height;ti3++)
      {
        for (int ti4=(ti3+1);ti4<width;ti4++)
        {
          tempinv1 = Recv2DInv[width*ti3+ti4];
          Recv2DInv[width*ti3+ti4] = Recv2DInv[width*ti3+ti4+(ti4-ti3)*(width-1)];
          Recv2DInv[width*ti3+ti4+(ti4-ti3)*(width-1)] = tempinv1;
        }
      }
    }

    Complex* Send1DInv = new Complex[width*RowCPU];
    if (rank==0)
    {
      for (int indi1=0;indi1<RowCPU*width;indi1++)
      {
        Send1DInv[indi1] = Recv2DInv[indi1];
      }
      for (int cpui2=1;cpui2<nCPUs;cpui2++)
      {
        rc = MPI_Send(&Recv2DInv[cpui2*width*RowCPU],RowCPU*width*sizeof(Complex),MPI_CHAR,cpui2,0,MPI_COMM_WORLD);
  // MPI send failure
        if (rc != MPI_SUCCESS)
        {
          cout << "Rank " << rank
               << " send failed, rc " << rc << endl;
          MPI_Finalize();
        }
      }
    }
    else
    {
      MPI_Status status;
      rc = MPI_Recv(&Send1DInv[0],RowCPU*width*sizeof(Complex),MPI_CHAR,0,0,MPI_COMM_WORLD, &status);
      if (rc != MPI_SUCCESS)
      {
        cout << "Rank " << rank
             << " recv failed, rc " << rc << endl;
        MPI_Finalize();
      }
    }

    int CurrentIndexInv1;
    Complex* Res1DInv = new Complex[width*RowCPU];
    for (int ii1=0;ii1<RowCPU;ii1++)
    {
      CurrentIndexInv1 = width*ii1; 
      Complex* Hi1 = new Complex[width];
  // extract the transposed columns that needs to perform fft
      for (int ji1=0;ji1<width;ji1++)
      {
        Hi1[ji1] = Send1DInv[CurrentIndexInv1+ji1];
      }
      Complex* hi1 = new Complex[width];
      Transform1D(Hi1,width,hi1,1);
  // save the outputs into a matrix
      for (int ki1=0;ki1<width;ki1++)
      {
        Res1DInv[CurrentIndexInv1+ki1] = hi1[ki1];
      }
    }

    Complex* Recv1DInv = new Complex[width*height];
  // other cpus send to cpu 0.
    if (rank!=0)
    {
      rc = MPI_Send(&Res1DInv[0],RowCPU*width*sizeof(Complex),MPI_CHAR,0,0,MPI_COMM_WORLD);
  // MPI send failure
      if (rc != MPI_SUCCESS)
      {
        cout << "Rank " << rank
             << " send failed, rc " << rc << endl;
        MPI_Finalize();
      }
    }
    else
    {
  // save cpu 0's own data.
      for (int ri1=0;ri1<RowCPU*width;ri1++)
      {
        Recv1DInv[ri1] = Res1DInv[ri1];
      }
  // cpu 0 receives from other cpus.
      int cpui2;
      for (cpui2=1;cpui2<nCPUs;cpui2++)
      {
        MPI_Status status;
        rc = MPI_Recv(&Recv1DInv[cpui2*width*RowCPU],RowCPU*width*sizeof(Complex),MPI_CHAR,cpui2,0,MPI_COMM_WORLD, &status);
        if (rc != MPI_SUCCESS)
        {
          cout << "Rank " << rank
               << " recv failed, rc " << rc << endl;
          MPI_Finalize();
        }
      }
      image.SaveImageData("MyAfterInverse.txt", Recv1DInv, width, height);
    }
}


// include DFT and inverse DFT
void Transform1D(Complex* h, int w, Complex* H, int flag)
{
  // Implement a simple 1-d DFT using the double summation equation
  // given in the assignment handout.  h is the time-domain input
  // data, w is the width (N), and H is the output array.
  if (flag==0)
  {
    for(int n=0;n<w;n++)
    {
      H[n] = 0.0;
      for(int k=0;k<w;k++)
      {
        Complex W(cos(2*M_PI*n*k/w),-sin(2*M_PI*n*k/w));
        H[n] = H[n] + W*h[k];
      }
    }
  }
  else
  {
    for (int i=0;i<w;i++)
    {
      H[i] = 0.0;
      for (int j=0;j<w;j++)
      {
        Complex W(cos(2*M_PI*i*j/w)/w,sin(2*M_PI*i*j/w)/w);
        H[i] = H[i]+W*h[j];
      }
    }
  }
}

int main(int argc, char** argv)
{
  string fn("Tower.txt"); // default file name
  if (argc > 1) fn = string(argv[1]);  // if name specified on cmd line
  // MPI initialization here
  int rc; 
  rc = MPI_Init(&argc,&argv);
  if (rc != MPI_SUCCESS) {
      printf ("Error starting MPI program. Terminating.\n");
      MPI_Abort(MPI_COMM_WORLD, rc);
      }
  Transform2D(fn.c_str()); // Perform the transform.
  //char *inputFN = "Tower.txt";
  MPI_Finalize();  // Finalize MPI here
}  
 
