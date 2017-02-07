PROG = SPFP
CC = g++ -std=c++11
CPPFLAGS = -g -O0 -Wall -fpermissive
LDFLAGS = -I/usr/local/include
LFLAGS = -L/usr/local/lib
LIBS = -lboost_thread -lboost_filesystem -lboost_program_options -lboost_system -lboost_date_time
OBJS = Group.o Utilities.o Cell.o Point.o Grid.o Value.o User.o Pair.o SPOs.o GPOs.o BasicGSQueries.o ArrayOperations.o CalculateProbability.o Entropy.o RenyiEntropy.o KatzScore.o WeightedEntropy.o SPFP.o

$(PROG) : $(OBJS)
	$(CC) $(LDFLAGS) -o $(PROG) $(OBJS) $(LFLAGS) $(LIBS)

-include $(OBJS:.o=.d)

SPFP.o : SPFP.cpp headersMemory.h
	$(CC) $(CPPFLAGS) -c SPFP.cpp
	$(CC) -MM SPFP.cpp > SPFP.d
Group.o : utilities/Group.cpp headers.h
	$(CC) $(CPPFLAGS) -c utilities/Group.cpp
	$(CC) -MM -c utilities/Group.cpp > Group.d
Utilities.o : utilities/Utilities.cpp headers.h
	$(CC) $(CPPFLAGS) -c utilities/Utilities.cpp
	$(CC) -MM utilities/Utilities.cpp > Utilities.d
Cell.o : GPOs/MemoryGrid/grid/Cell.cpp headersMemory.h
	$(CC) $(CPPFLAGS) -c GPOs/MemoryGrid/grid/Cell.cpp
	$(CC) -MM GPOs/MemoryGrid/grid/Cell.cpp > Cell.d
Grid.o : GPOs/MemoryGrid/grid/Grid.cpp headersMemory.h
	$(CC) $(CPPFLAGS) -c GPOs/MemoryGrid/grid/Grid.cpp
	$(CC) -MM GPOs/MemoryGrid/grid/Grid.cpp > Grid.d
Point.o : GPOs/MemoryGrid/grid/Point.cpp headersMemory.h
	$(CC) $(CPPFLAGS) -c GPOs/MemoryGrid/grid/Point.cpp
	$(CC) -MM GPOs/MemoryGrid/grid/Point.cpp > Point.d
SPOs.o : SPOs/MemoryMap/SPOs.cpp headersMemory.h
	$(CC) $(CPPFLAGS) -c SPOs/MemoryMap/SPOs.cpp
	$(CC) -MM SPOs/MemoryMap/SPOs.cpp > SPOs.d
GPOs.o : GPOs/MemoryGrid/GPOs.cpp headersMemory.h
	$(CC) $(CPPFLAGS) -c GPOs/MemoryGrid/GPOs.cpp
	$(CC) -MM GPOs/MemoryGrid/GPOs.cpp > GPOs.d
Value.o : SPOs/MemoryMap/Value.cpp  headersMemory.h
	$(CC) $(CPPFLAGS) -c SPOs/MemoryMap/Value.cpp
	$(CC) -MM SPOs/MemoryMap/Value.cpp > Value.d
User.o : SPOs/MemoryMapWeighted/User.cpp  headersMemory.h
	$(CC) $(CPPFLAGS) -c SPOs/MemoryMapWeighted/User.cpp
	$(CC) -MM SPOs/MemoryMapWeighted/User.cpp > User.d
Pair.o : SPOs/MemoryMapWeighted/Pair.cpp  headersMemory.h
	$(CC) $(CPPFLAGS) -c SPOs/MemoryMapWeighted/Pair.cpp
	$(CC) -MM SPOs/MemoryMapWeighted/Pair.cpp > Pair.d
BasicGSQueries.o : basicGSQueries/BasicGSQueries.cpp headers.h
	$(CC) $(CPPFLAGS) -c basicGSQueries/BasicGSQueries.cpp
	$(CC) -MM basicGSQueries/BasicGSQueries.cpp > BasicGSQueries.d
ArrayOperations.o : pTools/ArrayOperations.c pTools/ArrayOperations.h
	$(CC) $(cFLAGS) -c pTools/ArrayOperations.c
	$(CC) -MM pTools/ArrayOperations.c  > ArrayOperations.d
CalculateProbability.o : pTools/CalculateProbability.c pTools/CalculateProbability.h
	$(CC) $(cFLAGS) -c pTools/CalculateProbability.c
	$(CC) -MM pTools/CalculateProbability.c  > CalculateProbability.d
Entropy.o : pTools/Entropy.c pTools/Entropy.h
	$(CC) $(cFLAGS) -c pTools/Entropy.c
	$(CC) -MM pTools/Entropy.c  > Entropy.d
RenyiEntropy.o : pTools/RenyiEntropy.c pTools/RenyiEntropy.h
	$(CC) $(cFLAGS) -c pTools/RenyiEntropy.c
	$(CC) -MM pTools/RenyiEntropy.c  > RenyiEntropy.d
KatzScore.o : pTools/KatzScore.cpp headers.h
	$(CC) $(cFLAGS) -c pTools/KatzScore.cpp
	$(CC) -MM pTools/KatzScore.cpp  > KatzScore.d
WeightedEntropy.o : pTools/WeightedEntropy.c pTools/WeightedEntropy.h
	$(CC) $(cFLAGS) -c pTools/WeightedEntropy.c
	$(CC) -MM pTools/WeightedEntropy.c  > WeightedEntropy.d



.PHONY : clean
clean:
	rm -f core $(OBJS) *.d

.PHONY: cleanest
cleanest:
	rm -f core $(PROG) $(OBJS) *.d
