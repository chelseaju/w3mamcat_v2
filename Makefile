LIBS=-lcppcms -lbooster
CXX=g++

all: clean w3mamcat 

w3mamcat: app/w3mamcat.cpp w3mamcat_skin.cpp app/master.h model/*.cpp
	$(CXX) -shared -fPIC w3mamcat_skin.cpp model/*.cpp -o libw3mamcat_skin.dylib $(LIBS)
	$(CXX) -rdynamic app/w3mamcat.cpp model/*.cpp -o w3mamcat $(LIBS) 

w3mamcat_skin.cpp: view/master.tmpl view/contact.tmpl view/intro.tmpl view/input.tmpl view/result.tmpl app/master.h
	cppcms_tmpl_cc view/master.tmpl view/contact.tmpl view/intro.tmpl view/input.tmpl view/result.tmpl -o w3mamcat_skin.cpp

clean:
	rm -fr *.exe *.so w3mamcat_skin.cpp cppcms_rundir w3mamcat

