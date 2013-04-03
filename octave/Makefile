all:	support.oct

support.oct:	support.cc
	mkoctfile -I../support/Random123/include -o support.oct support.cc

.PHONY:	clean
clean:
	rm -f support.o support.oct

