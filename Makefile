all: synth decode

%.o: %.cpp
	$(CXX) -MMD -MP $(CPPFLAGS) $(CXXFLAGS) -o $@ -c $<

OBJS=decode.o synth.o synth_main.o

DEPS=$(OBJS:.o=.d)
-include $(DEPS)

decode: decode.o
	$(CXX) -o $@ $^ $(LDFLAGS)

synth: synth.o synth_main.o
	$(CXX) -o $@ $^ $(LDFLAGS)

clean:
	$(RM) synth $(OBJS) $(DEPS)
