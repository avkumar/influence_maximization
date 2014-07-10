CXXFLAGS =	-O2 -g -Wall -fmessage-length=0

OBJS =		Dijkstra.o mcMethods.o exactMethods.o GraphCPP.o

LIBS =

TARGET =  output.txt	

$(TARGET):	$(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS)

all:	$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)
