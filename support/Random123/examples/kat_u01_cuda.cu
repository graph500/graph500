/*
Copyright 2010-2011, D. E. Shaw Research.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

* Redistributions of source code must retain the above copyright
  notice, this list of conditions, and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions, and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

* Neither the name of D. E. Shaw Research nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#include "kat_u01_main.h"
#include "util_cuda.h"

#define KAT_KERNEL __global__
#define KAT_THREADID (blockDim.x * blockIdx.x + threadIdx.x)

#include "kat_u01_dev_execute.h"

void host_execute_tests(uint64_t *tests, size_t ntests, KatU01Result *results){
    CUDAInfo *infop;
    KatU01Result *results_dev;
    uint64_t *tests_dev;
    size_t tests_sz, results_sz;

    infop = cuda_init(NULL);
    tests_sz = ntests * sizeof(*tests_dev);
    CHECKCALL(cudaMalloc(&tests_dev, tests_sz));
    CHECKCALL(cudaMemcpy(tests_dev, tests, tests_sz, cudaMemcpyHostToDevice));
    
    results_sz = ntests * sizeof(*results);
    CHECKCALL(cudaMalloc(&results_dev, results_sz));

    printf("starting on %lu blocks with 1 threads/block\n", (unsigned long)ntests);
    fflush(stdout);

    dev_execute_tests<<<ntests, 1>>>(tests_dev, results_dev);

    CHECKCALL(cudaThreadSynchronize());
    CHECKCALL(cudaMemcpy(results, results_dev, results_sz, cudaMemcpyDeviceToHost));
    CHECKCALL(cudaFree(tests_dev));
    CHECKCALL(cudaFree(results_dev));
    cuda_done(infop);
}

