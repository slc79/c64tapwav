CXXFLAGS=-O2 -ffast-math -g -Wall

all: synth decode sync

%.o: %.cpp
	$(CXX) -MMD -MP $(CPPFLAGS) $(CXXFLAGS) -o $@ -c $<

OBJS=decode.o synth.o synth_main.o

DEPS=$(OBJS:.o=.d)
-include $(DEPS)

decode: decode.o
	$(CXX) -o $@ $^ $(LDFLAGS)

synth: synth.o synth_main.o
	$(CXX) -o $@ $^ $(LDFLAGS)

sync: sync.o
	$(CXX) -o $@ $^ $(LDFLAGS)

clean:
	$(RM) synth decode sync $(OBJS) $(DEPS)
