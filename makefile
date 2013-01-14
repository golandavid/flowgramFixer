#COPT =  -Wall -Wno-sign-compare -O3 -fPIC -fexceptions 
COPT =  -g -Wall -Wno-sign-compare -fPIC -fexceptions 



CC    = g++
CFLAGS = $(COPT) 


#.SUFFIXES:       # remove inbuilt definitions
#.SUFFIXES: .cpp .o # define default compilation procedure
#.cpp.o:            # .o files depend on .c files
#	$(CC) $(CFLAGS) $*.c # start with a tab character!!

# Default target
all: flowgramfixer
#	@echo "compilation done"

flowgramfixer: flowgramfixer.o
	$(CC) flowgramfixer.o -o flowgramfixer
flowgramfixer.o : flowgramfixer.cpp 
	$(CC) -c $(COPT) flowgramfixer.cpp -o flowgramfixer.o

# Target deleting unwanted files
clean:
	rm -f *.o *~ ut core 
