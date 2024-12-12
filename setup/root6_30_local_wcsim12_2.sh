#!/bin/bash

# This file contains the minimum configuration 
# to execute wcsim_root_to.C binary

module load DataManagement/xrootd/5.6.9
module load Programming_Languages/python/3.9.1
module load Analysis/root/6.30.06


wcsim_bin_path=/sps/t2k/eleblevec/WCSim-home/WCSim_FOR_READING_nuPRISM_first_HK_Hybrid_WITH_WCSIMROOT/WCSim_install/bin
source $wcsim_bin_path/this_wcsim.sh

export convert_bin=/sps/t2k/eleblevec/SimulationPackage/nuPRISM_hk_firstHybrid_WCSimRootToRoot/bin/new_wcsimroot_to_root


echo -e "\nWCSim binary path set to: $wcsim_bin_path"
echo -e "convert_bin env var set to :  $convert_bin\n"

# ------- REMARKS -------
# --- You should not modify this file



# -- Archive avant le 03/10/2024 --

# --- Le "module add Analysis/root/6.24.06" devrait aller chercher les trois librairies en dessous automatiquement 
# Programming_Languages/python/3.9.1
# Compilers/gcc/9.3.1
# DataManagement/xrootd/4.8.1
# --- Pour vérifier vous pouvez exécuter "module list" après avoir source ce .sh


# --- Le source exécute les 5 lignes suivantes :

# export WCSIM_BUILD_DIR=/sps/t2k/bquilain/HK/Reconstruction/official_2023/Gonzalo/WCSim-home/WCSim_build/mydir/
# export ROOT_INCLUDE_PATH=/sps/t2k/bquilain/HK/Reconstruction/official_2023/Gonzalo/WCSim-home/WCSim_build/mydir/include/WCSim:$ROOT_INCLUDE_PATH
# export PATH=/sps/t2k/bquilain/HK/Reconstruction/official_2023/Gonzalo/WCSim-home/WCSim_build/mydir/bin:$PATH
# export LD_LIBRARY_PATH=/sps/t2k/bquilain/HK/Reconstruction/official_2023/Gonzalo/WCSim-home/WCSim_build/mydir/lib:$LD_LIBRARY_PATH
# export WCSIMDIR=/sps/t2k/bquilain/HK/Reconstruction/official_2023/Gonzalo/WCSim-home/WCSim/