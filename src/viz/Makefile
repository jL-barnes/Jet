
EXE  = vjet

GLUT = /Library/Develop/CommandLineTools/SDKs/MacOSX.sdk/System/Library.Frameworks/GLUT.framework/Versions/A
HDF5 = /opt/local

LIBS = -L$(HDF5)/lib -framework GLUT -framework OpenGL -lXmu -lXext -lX11 -lm -lhdf5 -lz
INCL = -I$(GLUT)/Headers -I$(HDF5)/include

CC = gcc

default: $(EXE)

$(EXE): $(EXE).c
	$(CC) -c colorbar.c
	$(CC) -O3 -Wall -o $@ $< colorbar.o $(INCL) $(LIBS)

clean:
	rm -f *.o $(EXE)
