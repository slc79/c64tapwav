all: synth

%.o: %.cpp
	$(CXX) -MMD -MP $(CPPFLAGS) $(CXXFLAGS) -o $@ -c $<

OBJS=synth.o synth_main.o

DEPS=$(OBJS:.o=.d)
-include $(DEPS)

synth: synth.o synth_main.o
	$(CXX) -o $@ $^ $(LDFLAGS)

clean:
	$(RM) synth $(OBJS) $(DEPS)
