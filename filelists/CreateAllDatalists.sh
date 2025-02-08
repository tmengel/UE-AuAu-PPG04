#!/bin/bash
# rm -rf /sphenix/user/tmengel/UE-AuAu-PPG04/filelists/ana412_2023p015

# python3 CreateDataLists.py -b ana412 -d 2023p015 -r 23727
# python3 CreateDataLists.py -b ana412 -d 2023p015 -r 23735
# python3 CreateDataLists.py -b ana412 -d 2023p015 -r 23737
# python3 CreateDataLists.py -b ana412 -d 2023p015 -r 23738
# python3 CreateDataLists.py -b ana412 -d 2023p015 -r 23739
# python3 CreateDataLists.py -b ana412 -d 2023p015 -r 23740
# python3 CreateDataLists.py -b ana412 -d 2023p015 -r 23743
# python3 CreateDataLists.py -b ana412 -d 2023p015 -r 23745

cd /sphenix/user/tmengel/UE-AuAu-PPG04/filelists/ana412_2023p015
mkdir -p allruns
cp 237*/*.list allruns/.

# Get a list of all .list files
lists=$(ls allruns/*.list)

# Initialize an array for run numbers
runnumbers=()
output_file="data_job_args.list"
> "$output_file" # Clear the file if it already exists

for list in $lists; do
    runnumber=$(basename "$list" | sed -E 's/.*ana412_2023p015_000([0-9]+)-.*\.list/\1/')
    realpath=/sphenix/user/tmengel/UE-AuAu-PPG04/filelists/ana412_2023p015/allruns/$(basename "$list")
    echo "$runnumber, $realpath" >> "$output_file"
done

echo "Run numbers and paths written to $output_file"

# # Output the results
# echo "List files:"
# echo "$lists"
# echo "Run numbers:"
# echo "${runnumbers[@]}"
