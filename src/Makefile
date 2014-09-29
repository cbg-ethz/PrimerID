.PHONY: default

default: ../PID

OBJS = main.o reference.o alignment.o proper_read.o

../PID: $(OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ -o ../PID $(LIBS)

%.o : %.cpp $(wildcard %.hpp)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@
