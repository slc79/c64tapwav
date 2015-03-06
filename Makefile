CXXFLAGS=--std=gnu++0x -O2 -fno-math-errno -g -Wall
LDLIBS=-lavcodec -lavformat -lavutil -lswresample

all: synth decode sync level cleaner

%.o: %.cpp
	$(CXX) -MMD -MP $(CPPFLAGS) $(CXXFLAGS) -o $@ -c $<

OBJS=decode.o synth.o synth_main.o interpolate.o sync.o level.o

DEPS=$(OBJS:.o=.d)
-include $(DEPS)

decode: interpolate.o audioreader.o decode.o
	$(CXX) -o $@ $^ $(LDLIBS) $(LDFLAGS)

synth: synth.o synth_main.o
	$(CXX) -o $@ $^ $(LDFLAGS)

sync: interpolate.o sync.o
	$(CXX) -o $@ $^ $(LDFLAGS)

level: level.o
	$(CXX) -o $@ $^ $(LDFLAGS)

cleaner: cleaner.o
	$(CXX) -o $@ $^ $(LDFLAGS)

clean:
	$(RM) synth decode sync level cleaner $(OBJS) $(DEPS)
