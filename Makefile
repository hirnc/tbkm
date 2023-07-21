CXX = g++
CXXFLAGS = -O3 -fopenmp --std=c++11 -Wall -fpermissive -I. $(DEBUG)
LDLIBS += -lbamtools -lz -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5

all: create_read_count_matrix gff_coverage

clean:
	rm create_read_count_matrix gff_coverage *.o *.a

create_read_count_matrix:
	$(CXX) $(CXXFLAGS) $(INCLUDES) create_read_count_matrix.cpp gfflib.cpp hdf_base_depth_reader.cpp histd.cpp -o $@ $(LDLIBS)

gff_coverage:
	$(CXX) $(CXXFLAGS) $(INCLUDES) gff_coverage.cpp gfflib.cpp hdf_base_depth_reader.cpp histd.cpp -o $@ $(LDLIBS)

## dependency check ##
.KEEP_STATE:
.KEEP_STATE_FILE:.make.state.GNU-x86-Linux
