LIBS = -I/Users/laroque/Software/install/include -I/Users/laroque/Software/System_install/homebrew/include -L/Users/laroque/Software/install/lib -L/Users/laroque/Software/System_install/homebrew/lib -lgsl
#LIBS = -lgsl
ROOT = `root-config --cflags` `root-config --glibs`

aspen: tracker_aspen.cpp
	g++ -o aspen tracker_aspen.cpp $(LIBS) $(ROOT)
