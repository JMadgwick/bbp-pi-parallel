## bbp-pi-parallel - Calculate Pi Hex Digits in Parallel
This repo contains a simple program to calculate Pi Hex Digits in parallel, it uses the Bailey Borwein Plouffe formula.

It is a multi-threaded program which runs on CPUs and GPUs. It's an updated version of [bbpi](https://github.com/JMadgwick/bbpi/) which had not been updated in a while.

### How does it work?
Take a look at my [blog](http://madgwick.xyz/bbp-pi-parallel-calculation.php) where there's further information and a link to the paper containing the algorithms I used.

### What does it do?
By default it computes the 10 Millionth(10^7) hexadecimal digit of Pi (plus a few after that). And it doesn't take very long to do it either (see table below).

#### Usage
The CPU program optionally accepts two arguments, digit to calculate and number of threads to use (default all available). The GPU program accepts only digit to calculate.

#### Limitations
Both CPU and GPU versions of the program are limited by precision to calculating only the first 10^7 digits. Double precision (64-bit) floating point is used.
Most GPUs have very poor performance for 64-bit floating point operations compared with 32-bit operations. As such, the CPU program will usually be quicker. See my blog for more info.
For the CPU program, changing the code from `double` to `long double` allows 80-bit precision to be used (on x86). This will be slower than 64-bit, on AMD Zen 2 it runs 3 times slower.
Using 80-bit precision allows for the Pi Hex digits at 10^8 (one hundred million) to be calculated.

Most hardware doesn't support sufficient precision for digits greater than 10^8. Software can be used to implement greater precision (e.g. 128-bit), but this is very slow.
[GCC libquadmath](https://gcc.gnu.org/onlinedocs/libquadmath/) is one library which can be used for this, [GNU MPFR](https://www.mpfr.org/) is another.
A system which calculated 10^8 digits in 25 seconds (using `long double`), took an hour and a quarter to calculate the one billionth (10^9) digit (using `__float128`), which is 178x longer!

### Building

#### CPU
A compiler supporting C++11 is required. Use `march=native` to improve performance. I've also found clang produces quicker code than GCC.

`c++ -pthread -lm -std=c++11 -march=native -Ofast bbp-pi-parallel-cpu.cpp -o cpubbp.out`

#### GPU
The {HIPCC Compiler](https://github.com/ROCm-Developer-Tools/HIPCC) is required. See the [AMD ROCm Compiler Reference Guide](https://docs.amd.com/bundle/ROCm-Compiler-Reference-Guide-v5.5/page/Introduction_to_Compiler_Reference_Guide.html) for more information.
A supported GPU and its runtime is also required. For AMD this will be the HIP runtime, for Nvidia the propriatary driver and CUDA must be installed.

`hipcc bbp-pi-parallel-gpu.cpp -o gpubbp.out`

### The Files
+ bbp-pi-parallel-cpu.cpp
This is the code for the multi-threaded CPU implementation.
+ bbp-pi-parallel-gpu.cpp
This is the code for the HIP GPU implementation.

### Performance

These are some stats for CPUs/GPUs I have tested the program with.

| Digit | Time | Device |
| ---- | ------ | --------------- |
| 10^7 | 0.619 | Radeon VII |
| 10^7 | 0.587 | Ryzen 3700X |
| 10^7 | 2.444 | Xeon E5-2640 v3 |
| 10^8 | 14.542 | Ryzen 3700X |
| 10^8 | 23.214 | Xeon E5-2640 v3 |

The GPU program can (and should) be tuned for the characteristics for the GPU it's run on. This can be done by adjusting the `blocks`, `threadsPerBlock` & `perThreadRuns` variables.
