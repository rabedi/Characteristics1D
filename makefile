TOP_DIR= .
SHELL=/bin/sh
export ARCH = KL_config
include $(TOP_DIR)/$(ARCH)
SRCFILES := $(wildcard src/*.cpp src/*.cc)

.PHONY: main 

all: main

main:  
	$(CXX) $(CXX_FLAGS) $(CFLAGS)  $(INCLUDES)  main.cpp   -o solver  $(LIBS) $(SRCFILES)



clean:
	( rm main )
	
