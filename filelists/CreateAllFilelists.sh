#!/bin/bash
echo "Creating data lists for run 23745"
python3 CreateDataLists.py -b ana410 -d 2023p014 -r 23745

echo "Create jet30 lists"
mkdir -p sim_run19_type11_jet30_nop 
cd sim_run19_type11_jet30_nop
CreateFileList.pl -type 11 -run 19 -nop DST_CALO_CLUSTER DST_GLOBAL DST_MBD_EPD DST_TRUTH_JET

echo "Create jet10 lists"
cd ..
mkdir -p sim_run19_type12_jet10_nop
cd sim_run19_type12_jet10_nop
CreateFileList.pl -type 12 -run 19 -nop DST_CALO_CLUSTER DST_GLOBAL DST_MBD_EPD DST_TRUTH_JET

echo "Create hijing MB lists"
cd ..
mkdir -p sim_run19_type4_hijing_nop
cd sim_run19_type4_hijing_nop
CreateFileList.pl -type 4 -run 19 -nop DST_CALO_CLUSTER DST_GLOBAL DST_MBD_EPD