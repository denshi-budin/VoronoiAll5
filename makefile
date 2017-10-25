IDIR =include
CC=gcc
CFLAGS=-I$(IDIR) -std=gnu99

SDIR=src

ODIR=obj
LDIR =lib

LIBS=-lm

_DEPS = atom.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = MDmain.o Placement.o getrandom.o Sort.o Euler.o InOrOut.o writexyz.o Timer.o inputpoint.o 
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: $(SDIR)/%.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

voronoi: $(OBJ)
	gcc -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ voronoi*