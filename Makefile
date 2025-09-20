CXXFLAGS = -std=c++11 -O3

sunsky: SunSky.cpp SunSky.hpp SunSkyTool.cpp
	$(CXX) $(CXXFLAGS) -o $@ SunSky.cpp SunSkyTool.cpp

clean:
	$(RM) sunsky
