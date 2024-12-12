## Requirements :

- a compiled wcsim
- the root that was used to compile wcsim

## Usage :

1. Modify **root6_30_local_wcsim12_2.sh** with you path to your wcsim (the wcsim_bin_path variable) then source it

```bash
source setup/root6_30_local_wcsim12_2.sh
```

1. Compile the code 

```bash
make
```

1. Define **debug_config.txt** as you need & execute wcimsroot_to_root

```bash
bin/new_wcsimroot_to_root debug_config.txt
```

# Potential errors

⚠️ You might encounter an error with the environment variable **WCSIMDR**

It is due to the fact that the wcsim developpers changed **WCSIMDIR** to **WCSIM_SOURCE_DIR** in their **this_wcsim**.**sh**

In this case you have to modify the MakeFile of wcsimroot_to_root this way : 

```bash
# BEFORE
# ------- WCSIM -------- #
WC_INC = $(WCSIMDIR)/include
WC_SRC= $(WCSIMDIR)/src
CXXFLAGS += -I$(WC_SRC) -I$(WC_INC) 

# AFTER
# ------- WCSIM -------- #
WC_INC = $(WCSIM_SOURCE_DIR)/include
WC_SRC= $(WCSIM_SOURCE_DIR)/src
CXXFLAGS += -I$(WC_SRC) -I$(WC_INC) 
```

Then you should be able to run (after source setup/…)

```bash
make
```

to compile wcsimroot_to_root


## GetVtx_old(k)

Needs to be change to GetVtx(k)