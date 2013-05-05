CXXFLAGS=-O2 -ffast-math -g -Wall

all: synth decode sync level

%.o: %.cpp
	$(CXX) -MMD -MP $(CPPFLAGS) $(CXXFLAGS) -o $@ -c $<

OBJS=decode.o synth.o synth_main.o interpolate.o sync.o level.o

DEPS=$(OBJS:.o=.d)
-include $(DEPS)

decode: interpolate.o decode.o
	$(CXX) -o $@ $^ $(LDFLAGS)

synth: synth.o synth_main.o
	$(CXX) -o $@ $^ $(LDFLAGS)

sync: interpolate.o sync.o
	$(CXX) -o $@ $^ $(LDFLAGS)

level: level.o
	$(CXX) -o $@ $^ $(LDFLAGS)

clean:
	$(RM) synth decode sync level $(OBJS) $(DEPS)
