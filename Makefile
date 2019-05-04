include Makefile.in

OPTFLAGS = -O2

ifneq ($(NETCDF4HOME),)
  NETCDFLIBDIR=-L$(NETCDF4HOME)/lib
  NETCDFINCLUDE=-I$(NETCDF4HOME)/include
  NETCDFLD=-lnetcdf
  NETCDFDEF = -DUSENETCDF
  NETCDFSRC= mynetcdf.c
else
  NETCDFLIBDIR = 
  NETCDFINCLUDE = 
  NETCDFLD = 
  NETCDFDEF = -DNONETCDF
  NETCDFSRC= mynetcdf-nonetcdf.c
endif

#LIBS = $(NETCDFLD)
LIBDIR = $(NETCDFLIBDIR)
#LDFLAGS = -lm $(LIBDIR) $(LIBS)
INCLUDES = $(NETCDFINCLUDE) 
DEFINES = $(NETCDFDEF)
CFLAGS = $(INCLUDES) #$(DEFINES)

LIBS = -lnetcdf -lm #-o pgm

CC = gcc

all: transport_model.c
	${CC} -o transport_model transport_model.c functions.c $(CFLAGS) $(LIBDIR) $(LIBS)

clean:
	rm -f *.o

transport_model.o: functions.h
