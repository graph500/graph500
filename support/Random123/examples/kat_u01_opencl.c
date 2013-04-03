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

#include "util_opencl.h"
#include <kat_u01_opencl_kernel.i>

void host_execute_tests(uint64_t *tests, size_t ntests, KatU01Result *results)
{
    UCLInfo *infop;
    cl_kernel kern;
    size_t tests_sz, results_sz;
    cl_mem tests_dev, results_dev;
    const char *kernelname = "dev_execute_tests";
    cl_int err;

    infop = opencl_init(NULL, opencl_src, "");
    CHECKERR(kern = clCreateKernel(infop->prog, kernelname, &err));
    infop->wgsize = 1;

    tests_sz = ntests * sizeof(*tests);
    CHECKERR(tests_dev = clCreateBuffer(infop->ctx, CL_MEM_READ_WRITE|CL_MEM_USE_HOST_PTR, tests_sz, tests, &err));
    CHECK(clEnqueueWriteBuffer(infop->cmdq, tests_dev, CL_TRUE, 0, tests_sz, tests, 0, 0, 0));
    
    results_sz = ntests * sizeof(*results);
    CHECKERR(results_dev = clCreateBuffer(infop->ctx, CL_MEM_WRITE_ONLY, results_sz, 0, &err));
    
    CHECK(clSetKernelArg(kern, 0, sizeof(cl_mem), (void*)&tests_dev));
    CHECK(clSetKernelArg(kern, 1, sizeof(cl_mem), (void*)&results_dev));
    printf("queuing kernel for %lu threads with %lu work group size\n",
	   (unsigned long)ntests, (unsigned long)infop->wgsize);
    fflush(stdout);
    CHECK(clEnqueueNDRangeKernel(infop->cmdq, kern, 1, 0, &ntests, &infop->wgsize, 0, 0, 0));
    CHECK(clFinish(infop->cmdq));

    CHECK(clEnqueueReadBuffer(infop->cmdq, results_dev, CL_TRUE, 0, results_sz, results, 0, 0, 0));

    CHECK(clReleaseMemObject(tests_dev));
    CHECK(clReleaseMemObject(results_dev));
    CHECK(clReleaseKernel(kern));
    opencl_done(infop);
}
