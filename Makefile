CPP 		= g++
CXXFLAGS	= -g -O3 -Wall -fPIC -D_REENTRANT -Wno-deprecated -fpermissive -std=c++17
LIBS=

# ------- ROOT --------- #
ROOT_CXXFLAGS	:= $(shell root-config --cflags)
CXXFLAGS        += $(ROOT_CXXFLAGS)

ROOT_LIBS    := $(shell root-config --libs)
ROOT_GLIBS   := $(shell root-config --glibs)
LIBS +=  $(ROOT_LIBS) $(ROOT_GLIBS)

# ------- WCSIM -------- #
WC_INC = $(WCSIMDIR)/include
WC_SRC= $(WCSIMDIR)/src
CXXFLAGS += -I$(WC_SRC) -I$(WC_INC) 

WCSIM_LIBS := -L$(WCSIM_BUILD_DIR)/lib -lWCSimRoot
LIBS += $(WCSIM_LIBS)

# ------- DISPLAY FLAGS -------- # 06/01/2024 - Erwan : Pour l'instant pas utile, mais à voir
# ifeq ($(VERBOSE),1)
# CXXFLAGS += -DVERBOSE
# endif


# ------- COMPILATION LOOP -------- #
TARGET      = new_wcsimroot_to_root
SRC_DIR     = src
OBJ_DIR     = obj
OUTPUT_DIR  = bin

SRC_FILES   = $(wildcard $(SRC_DIR)/*.cpp)
OBJ_FILES   = $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRC_FILES))

all: $(OUTPUT_DIR)/$(TARGET)

$(OUTPUT_DIR)/$(TARGET): $(OBJ_FILES) | $(OUTPUT_DIR)
	@echo "Now make $@"
	@echo $(LIBS)
	
	@$(CPP) -o $@ $^ $(CXXFLAGS) $(LIBS)
	
	@echo "..Compile done! "

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	@echo "$<"
	@echo "Start Compiling $<"

	@$(CPP) $(CXXFLAGS) -c $< -o $@
	
	@echo ".. Compiling Object Files $<   --> done"
	@echo ""

$(OBJ_DIR):
	@mkdir -p $(OBJ_DIR)

$(OUTPUT_DIR):
	@mkdir -p $(OUTPUT_DIR)

clean: 
	@echo "Now Clean Up"
	rm -f $(OUTPUT_DIR)/$(TARGET) *~ $(SRC_DIR)/*~ $(OBJ_DIR)/*.o core

# $(OUTPUT_DIR)/$(TARGET) *~ cleans the binaries in the bin/
# $(SRC_DIR)/*~				 cleans the binaries in the src
# $(OBJ_DIR)/*.o			 cleans all the files finishing with .o in obj
# core : clean the rest of the useless files generated for the compilation purpose noky