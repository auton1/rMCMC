CPP = g++
EXECUTABLE = rMCMC

CFLAGS = -O2 -m64 
CPPFLAGS = -O2 -D_FILE_OFFSET_BITS=64 -std=c++0x
LIB = -lz 

OBJS = rMCMC.o LogFactorial.o event.o interval.o \
	readFile.o set_intervals.o read_nbmeioses.o\

rMCMC: $(OBJS)
	$(CPP) $(CPPFLAGS) $(OBJS) -o rMCMC $(LIB)

# pull in dependency info for *existing* .o files
-include $(OBJS:.o=.d)

%.o: %.cpp
	$(CPP) -c $(CPPFLAGS) $*.cpp -o $*.o
	$(CPP) -MM $(CPPFLAGS) $*.cpp > $*.d

# remove compilation products
clean:
	@rm -f rMCMC *.o *.d

