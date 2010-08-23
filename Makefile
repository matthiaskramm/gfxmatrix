CC=gcc -Wimplicit -O3  -fomit-frame-pointer
OPT=-mfpmath=sse -msse2 -mcpu=pentium4 -march=pentium4 -mtune=pentium4 -O3 -fomit-frame-pointer

OBJS=texture.o image.o filter.o filter_steerpyr.o gsltools.o stat.o matrix.o corr.o gradopt.o graphcut.o filter_misc.o video.o filter_wavelet.o segment.o kmeans.o txtsynth.o
LIBS=-lz -lgsl -lgslcblas -lm -lX11 -lXm -lpthread -lxml2 -lfftw3

all: color imglib.a test corr matrix show_filter graphcut testwindow analysetexture edge video gfximage.so

color: color.c
	$(CC) -DHAVE_GSL color.c -o color ../png.o -lgsl -lgslcblas  -lz

edge: edge.c imglib.a matrix.h
	$(CC) -DHAVE_GSL edge.c -o edge imglib.a gfxwindow.o $(LIBS)

%.o: %.c
	$(CC) -DHAVE_GSL -c $< -o $@

gfxwindow.o: gfxwindow.c gfxwindow.h gfxwindow_unix.c gfxwindow_win32.c Makefile
	gcc -c -g -O2 gfxwindow.c -o gfxwindow.o

filter.o: filter.c filter.h
matrix.o: matrix.c
	$(CC) -DHAVE_GSL -march=pentium3 -mtune=pentium3 -fno-strict-aliasing  -c $< -o $@
image.o: image.c
video.o: video.c video.h
gradopt.o: gradopt.c
graphcut.o: graphcut.c
stat.o: stat.c
gsltools.o: gsltools.c
corr.o: corr.c corr.h
texture.o: texture.c texture.h stat.h filter.h
filter_steerpyr.o: filter_steerpyr.c filter.h
filter_misc.o: filter_misc.c filter.h
filter_wavelet.o: filter_wavelet.c filter.h
segment.o: segment.c segment.h

kmeans.o: kmeans.asm kmeans.h
	nasm -f elf kmeans.asm -o kmeans.o

video: testvideo.c imglib.a
	$(CC) testvideo.c -o video gfxwindow.o imglib.a $(LIBS)

testwindow: testwindow.c gfxwindow.o
	$(CC) testwindow.c -o $@ gfxwindow.o $(LIBS) ../png.o

test: test.c imglib.a
	$(CC) -DHAVE_GSL -o $@ test.c imglib.a $(LIBS)

corr: testcorr.c imglib.a
	$(CC) -DMAIN -DHAVE_GSL -o $@ testcorr.c imglib.a $(LIBS)

matrix: matrix.c image.o gsltools.o stat.o
	$(CC) $(OPT) -DHAVE_GSL -DMAIN matrix.c image.o stat.o gradopt.o gsltools.o -o matrix $(LIBS)

show_filter: show_filter.c imglib.a
	$(CC) -I /usr/include/libxml2/ show_filter.c -o show_filter imglib.a $(LIBS)

graphcut: graphcut.c graphcut.h
	gcc -O3 -fomit-frame-pointer -DMAIN -DHAVE_GSL graphcut.c -o graphcut imglib.a $(LIBS)

f: f.c
	gcc -O3 -fomit-frame-pointer -DMAIN -DHAVE_GSL f.c -o f imglib.a $(LIBS)

imglib.a: $(OBJS)
	ar cru imglib.a $(OBJS)

analysetexture: analysetexture.c image.h imglib.a
	$(CC) -DHAVE_GSL analysetexture.c imglib.a $(LIBS) -o analysetexture

txtsynth.o: txtsynth.c Makefile
	$(CC) -I../lib -c txtsynth.c -o txtsynth.o

txtsynth: txtsynth_main.c txtsynth.o ../lib/imglib.a  ../lib/gfxwindow.o
	$(CC) -DMAIN -I ../lib/ -Wimplicit txtsynth_main.c txtsynth.o -o txtsynth imglib.a gfxwindow.o -lm -lz -lgsl -lgslcblas -lfftw3 $(LIBS)

gfximage.o: gfximage.c
	$(CC) -fPIC -c $(PYTHON_INCLUDES) gfximage.c -o gfximage.o

gfximage.so: gfximage.o $(OBJS) Makefile txtsynth.o
	$(CC) $(PYTHON_INCLUDES) -shared $(OBJS) gfximage.o -o gfximage.so $(PYTHON_LIB) $(LIBS)

clean:
	rm -f rev*.dat test*.dat rev*.png test*.dat $(OBJS) imglib.a
