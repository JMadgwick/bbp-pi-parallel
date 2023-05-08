// Parallel Implementation of the Bailey–Borwein–Plouffe Formula for Pi (Improved Program)
// Copyright (C) 2023 J. Madgwick

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <iostream>
#include <cmath>
// Header files for the HIP API
#include <hip/hip_runtime.h>
#include <hip/hip_runtime_api.h>

// Left-Right Binary algorithm for exponentiation with Modulo 16^n mod k
__host__ __device__ double expoMod(double n, double k)
  {
    static int init = 0; // Store whether table initialised
    static int pwrtbl[32]; // Table to store powers of 2
    static int highestpwrtblpwr = 0; // After init, points to largest value in the table, which is <= n
    while (!init) // Find largest power needed, as n counts down, do this only once
    {
      pwrtbl[0] = 1; // In theory first can be 0 but in practice this function is only called if n > 0, in which case greatest power is 1 or more
      while (pwrtbl[highestpwrtblpwr] < n) // Used instead of for loop in order to check before iterating
      {
        pwrtbl[highestpwrtblpwr+1] = pwrtbl[highestpwrtblpwr] * 2; // Only generates up to and including n, nothing larger
        highestpwrtblpwr++;
      }
      init = 1; // Set table as initialised
    }

    int bitsneeded = 0; // The number of binary positions and the location in the table of the highest needed power
    for (; pwrtbl[bitsneeded]<n;bitsneeded++); // Move through table to find position of the largest power of two < n
    if (pwrtbl[bitsneeded]!=n){bitsneeded--;} // Because the increment is applied before the condition check, we need to decrement by one, unless the item was equal to n

    // First set t to be the largest power of two such that t ≤ n, and set r = 1.
    int t = pwrtbl[bitsneeded];
    double r = 1;

    // Loop by the number of Binary positions
    for (int i = 0; i <= bitsneeded; i++)
    {
      if (n>=t) // if n ≥ t then
      {
        r = r * 16.0; // r ← br
        r = r - static_cast<int>(r/k)*k; // r ← r mod k
        n = n - t; // n ← n − t
      }
      t = t/2; // t ← t/2
      if (t>=1) // if t ≥ 1 then
      {
        r = r * r; // r ← r^2
        r = r - static_cast<int>(r/k)*k; // r ← r mod k
      }
    }
    return r;
  }

// Left Portion for one thread - just does one term
__device__ double leftPortionThreaded(int terms, int k, int j, int d)
{
  int kinit = k;
  double s = 0;
  for (;k < (kinit+terms);k++)
  {
    double numerator,denominator;
    denominator = 8 * k + j;
    numerator = expoMod(d - k, denominator); // Binary algorithm for exponentiation with Modulo must be used becuase otherwise 16^(d-k) can be very large and overflows
    s = s + numerator/denominator;
    s = s - std::floor(s);
  }
  return s;
}

__global__ void kern(double *gpu_threadResults, int perThreadRuns, int j, int d, int k){
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int k2 = k+(idx*perThreadRuns); // Thread 0 starts at k+0, thread 1 at k+(1*perThreadRuns) &c.
  gpu_threadResults[idx] = leftPortionThreaded(perThreadRuns,k2,j,d);
}

// Bailey–Borwein–Plouffe Formula 16^d x Sj
double bbpf16jsd(int j, int d)
  {
    double s = .0;
    double numerator,denominator;
    double term;
    int blocks = 80;
    int threadsPerBlock = 60;
    int noOfThreads= blocks * threadsPerBlock;
    int perThreadRuns = 2000;
    double threadResults[noOfThreads]; // For storing results from threads
    double *gpu_threadResults;
    hipMalloc(&gpu_threadResults,noOfThreads*sizeof(double));

    // Left Portion
    int k = 0;
    while (k < d)
    {
      if (k + (perThreadRuns*noOfThreads) < d) // Only make threads for k up to less than d
      {
        // Launch HIP Kernel, then copy results from GPU to CPU
        hipLaunchKernelGGL(kern,dim3(blocks),dim3(threadsPerBlock),0,0,gpu_threadResults,perThreadRuns,j,d,k);
        hipMemcpy(threadResults,gpu_threadResults,noOfThreads*sizeof(double),hipMemcpyDefault);

        k = k + perThreadRuns*noOfThreads; // We need to run (perThreadRuns) results on each thread because the thread overhead is much to great to run just 1
        for (int i2 = 0; i2 < noOfThreads;i2++)
        {
          s = s + gpu_threadResults[i2];
          s = s - static_cast<int>(s);
        }
      } else // If we are almost done and k + noOfThreads > d then do the last few terms single threaded
      {
        denominator = 8 * k + j;
        numerator = expoMod(d - k, denominator); // Binary algorithm for exponentiation with Modulo must be used becuase otherwise 16^(d-k) can be very large and overflows
        term = numerator/denominator;
        s = s + term;
        s = s - static_cast<int>(s);
        k++;
      }
    }
    // Right Portion
    for (int k = d; k <= d+100; k++)
    {
      numerator = pow(16, d - k);
      denominator = 8 * k + j;
      term = numerator/denominator;
      if (term<1e-17) {break;}
      s = s + term;
      s = s - static_cast<int>(s);
    }
    return s;
  }

// Bailey–Borwein–Plouffe Formula Calculation
void bbpfCalc(double *pidec,int *place)
  {
    int tempn = *place;
    double result;
    result = (4.*bbpf16jsd(1,tempn))-(2.*bbpf16jsd(4,tempn))-bbpf16jsd(5,tempn)-bbpf16jsd(6,tempn);
    result = result - static_cast<int>(result) + 1.;
    *pidec = result;
  }

void toHex(char *out, double *in)
{
  char hexNumbers[] = "0123456789ABCDEF";

  for (int i = 0; i <= 8;i++)
  {
    *in = 16.0 * (*in - std::floor(*in));
    out[i] = hexNumbers[static_cast<int>(*in)];
  }
}

int main(int argc, char *argv[]) {
  std::cout << "Bailey–Borwein–Plouffe Formula for Pi" << std::endl;
  std::cout << "Built: " << __DATE__ << " " << __TIME__ << " with HIP Version: " << HIP_VERSION_MAJOR << "." << HIP_VERSION_MINOR << "." << HIP_VERSION_PATCH << std::endl << std::endl;
  hipDeviceProp_t GPUdevice; hipGetDeviceProperties(&GPUdevice, 0);
  std::wcout << "-------- Detected HIP Device Details --------" << std::endl
  << "          Name: " << GPUdevice.name << std::endl
  << "     Total RAM: " << GPUdevice.totalGlobalMem/pow(1024,2) << " (MB)" << std::endl // RAM is shown is MB output from API is bytes
  << " Compute Units: " << GPUdevice.multiProcessorCount << std::endl << std::endl;
  int placeNo = (argc >= 2) && (std::atoi(argv[1]) > 0) ? std::atoi(argv[1]) - 1 : 10000000 - 1; // Accurate to 10000000
  std::cout << "Calculating Position: " << (placeNo + 1) << std::endl;
  double piDec;
  bbpfCalc(&piDec, &placeNo);
  char hexOutput[] = "000000000";
  toHex(hexOutput, &piDec);
  std::cout << "Pi Estimation Hex: " << hexOutput << std::endl;
}
