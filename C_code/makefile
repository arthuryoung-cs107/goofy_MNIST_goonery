CC=gcc-9
OS=OSX
CFLAGS= -lm -Ofast -fomit-frame-pointer -ffast-math -fstrict-aliasing -finline-functions -I. -DGSL

LDIR= -I/usr/local/include/
LINK= -L/usr/local/lib/ -lgsl -lgslcblas

DEPS3 = nrutil.h aux_functions.h GSL_tomfoolery.h line_search.h

OBJS3 = MNIST_gaussian_denoiser.o nrutil.o aux_functions.o GSL_tomfoolery.o line_search.o

%.o: %.c $(DEPS3)
	$(CC) -c -o $@ $< $(CFLAGS) $(LDIR)

run_denoiser: $(OBJS3)
	$(CC) -o $@ $^ $(CFLAGS) $(LINK)


clean:
	rm -f  run_denoiser
	rm -f *.o

clean_dat:
	rm -f *.dat
clean_csv:
	rm -f ./MNIST_interim/*.csv
