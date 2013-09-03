no_target_specified: runcore
	@echo
	@echo The default make rule is equivalent to \'make runcore\' which runs only the most basic tests.
	@echo The following \'meta-targets\' are available:
	@echo "      " $(meta_targets)
	@echo Here is the complete list of individual program targets:
	@echo "      " $(all_primary_targets)
	@echo Prepend \'run\' to any of the program targets or metatargets
	@echo to run the binary and check for a zero exit status.
	@echo Adding force=1 on the command line causes all targets to be considered out-of-date.
.PHONY: no_target_specified

# metatargets are variables which get mapped by METATARGET_template
meta_targets:=kat core aesni timing c cpp gsl thread cuda opencl

# Platform metatargets: each one typically has specific requirements in the build environment.
# c is C99 (will work in MSVC), cpp is C++98, gsl requires the GNU Scientific Library
# (specifically, the gsl-config program in the PATH), thread requires POSIX threads,
# CUDA requires NVIDIA CUDA 3.x or newer, OpenCL requires OpenCL includes & libraries
# (e.g. AMD APP SDK, NVIDIA SDK)
c:=kat_c pi_capi simple ut_ars time_serial
cpp:=kat_cpp pi_uniform pi_cppapi simplepp ut_carray ut_M128 ut_features ut_ReinterpretCtr ut_Engine ut_aes pi_aes timers pi_microurng time_boxmuller ut_ua ut_uniform
gsl:=pi_gsl ut_gsl
thread:=time_thread
cuda:=pi_cuda pi_cudapp time_cuda kat_cuda time_boxmuller_cuda
opencl:=pi_opencl time_opencl kat_opencl

# Convenience metatargets: these are to help developers test functional subsets across platforms
kat:=kat_c kat_cpp
core:=$(c) $(cpp)
aesni:=pi_aes ut_aes ut_ars
timing:=timers time_serial time_thread

$(gsl) : override LDLIBS += `gsl-config --libs`
$(gsl) : override CFLAGS += `gsl-config --cflags`

$(opencl) : % : %_kernel.i
$(opencl) : override LDLIBS+=-lOpenCL
$(opencl) : override CFLAGS+=-I.
# Note, the Intel OpenCL SDK (1.5) has unresolved C++ symbols in its
# libOpenCL.so Even though 'main' is a C program, you may need to link
# it with a C++ compiler-driver, e.g., g++.  Since this Makefile does
# compile-and-link in one step, use something like:
# $(opencl) : CC=g++ -xc
# which will invoke the g++ compiler-driver, but will treat the
# program as C rather than C++.

# N.B. gcc -pthread (without the -l) on linux at compile time also
# adds -D_REENTRANT.  Unfortunately -pthread is unrecognized by
# SunPRO.  Posix says that -lpthread is portable (as if anyone cares
# what Posix says), and it seems to work in all the environments we've
# tested.  Gcc's features.h says that _THREAD_SAFE is "often used by
# other systems" as a synonym for _REENTRANT.  Cross your fingers...
$(thread) : override LDLIBS+=-lpthread
$(thread) : override CPPFLAGS+=-D_REENTRANT=1 -D_THREAD_SAFE=1

all_primary_targets += $(addsuffix _kernel.i, $(opencl))

################################################
# Generic boilerplate from here down:
vpath %.c $(srcdir/)
vpath %.cpp $(srcdir/)
vpath %.cu $(srcdir/)
vpath %.ocl $(srcdir/)

define METATARGET_template
.PHONY: $(1)
$(1) : $(filter-out $(SKIP_TARGETS), $($(1)))
.PHONY: run$(1)
run$(1) : $(addprefix run, $(filter-out $(SKIP_TARGETS), $($(1))))
all_primary_targets += $($(1))
endef

$(foreach T,$(meta_targets), $(eval $(call METATARGET_template,$(T))))

# sort also does 'uniq'
all_primary_targets:=$(sort $(all_primary_targets))

INC=$(srcdir/)../include
override CPPFLAGS += -I$(INC)

ifndef NVCC
NVCC:=nvcc
endif
# The rngs are *very* slow without optimization.  In the simplest case,
# where the user just calls 'make', we don't want them to see terrible
# performance.  Unfortunately, this might surprise someone
# who says, e.g., make CPPFLAGS=-O0.  Oh well...
ifndef CFLAGS
CFLAGS:=-O
endif
ifndef CXXFLAGS
CXXFLAGS:=-O
endif

%.i : %.ocl
	CC="$(CC)" CPPFLAGS="$(CPPFLAGS)" $(srcdir/)./gencl.sh $< $@

% : %.cu
	$(NVCC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) $(TARGET_ARCH) $< $(LOADLIBES) $(LDLIBS) -o $@

% : %.c
	$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) $(TARGET_ARCH) $< $(LOADLIBES) $(LDLIBS) -o $@

% : %.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) $(TARGET_ARCH) $< $(LOADLIBES) $(LDLIBS) -o $@

run% : %
	./$^ $(RUN_ARGS)

# In lieu of autodepends, just say that all the compilation targets depend on all the headers.
hdrs:=$(wildcard $(srcdir/)*.h $(srcdir/)gsl/*.c $(INC)/Random123/*.h $(INC)/Random123/*.hpp $(INC)/Random123/*/*.h $(INC)/Random123/*/*.hpp)
misc:=$(wildcard $(srcdir/)*.cu $(srcdir/)*.ocl)
$(all_primary_targets) : $(hdrs)
$(misc) : $(hdrs)

# If you put force=y on the command line, then $(all_primary_targets) will be
# depend on FORCE, and hence will not be up-to-date.
ifdef force
$(all_primary_targets) : FORCE
FORCE:
endif

.PHONY : echo_primary_targets
echo_primary_targets:
	@echo $(all_primary_targets)

.PHONY : clean veryclean
clean:
	rm -f $(all_primary_targets)

veryclean:
	rm -f $(all_primary_targets) *.o \#* *~ *.pdb *.exe *.obj *.ilk *.suo
