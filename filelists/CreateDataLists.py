#!/bin/python3
import os
import sys
import argparse

def config_output_directory(anabuild, database):

    # set up the directory structure
    current_dir = os.getcwd()
    output_dir = os.path.join(current_dir, f'{anabuild}_{database}')

    # clear the directory if it exists
    # if os.path.exists(output_dir):
        # os.system(f'rm -rf {output_dir}')
    os.makedirs(output_dir, exist_ok=True)
    
    return output_dir

def split_run_list(runlist, output_directory, n_file_per_list=10):

    # get the base name of the runlist file
    runlist_base = os.path.basename(runlist)
    runlist_base = runlist_base.split(".")[0]

    # make sure runlist file exists
    if not os.path.exists(runlist):
        print("runlist file does not exist: ", runlist)
        sys.exit(1)

    # read first line of the runlist file
    with open(runlist, "r") as f:
        first_line = f.readline().strip()
    
    # get the run number
    run_number = first_line.split("-")[1]
    run_number_for_dir = str(int(run_number))

    # get production name
    detector = runlist_base.split("-")[0]
    # make detector name all uppercase
    detector = detector.upper()


    data_type = first_line.split("_")[2]
    production_name = first_line.split("_")[3]
    year = first_line.split("_")[4]
    year = year.split("-")[0]

    directory_name = f'{run_number_for_dir}'
    directory_name = os.path.join(output_directory, directory_name)
    # clear the directory if it exists
    if os.path.exists(directory_name):
        os.system(f'rm -rf {directory_name}')
    
    os.makedirs(directory_name, exist_ok=True)

    filebase = f'{detector}_{data_type}_{production_name}_{year}_{run_number}-'
    filepath = os.path.join(directory_name, filebase)
    
    # create new runlist based on the number of files per list
    current_nfile = -1
    
    # read all lines in the runlist file
    with open(runlist, "r") as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            if i % n_file_per_list == 0:
                if current_nfile > 0:
                    new_runlist.close()
                current_nfile += 1
                # new file is filepath{%04d}
                new_runlist = open(f'{filepath}{current_nfile:04d}.list', "w")
            new_runlist.write(line)
        new_runlist.close()
    
def gen_run_list(anabuild, database, runnumber=-1, n_file_per_list=10):
    
    # generate the runlist file
    output_dir = config_output_directory(anabuild, database)

    # cd to the output directory
    os.chdir(output_dir)

    # create output log file
    log_file_path = os.path.join(output_dir, f"{anabuild}_{database}_runlist_gen.log")

    # create directory for the runlists
    runlist_dir = os.path.join(output_dir, "runlist_tmp")
    os.makedirs(runlist_dir, exist_ok=True)

    # cd to the runlist directory
    os.chdir(runlist_dir)

    with open(log_file_path, "w") as f:
        f.write(f"anabuild: {anabuild}\n")
        f.write(f"database: {database}\n")
        f.write(f"runnumber: {runnumber}\n")
        f.write(f"output_dir: {output_dir}\n")
        f.write(f"runlist_dir: {runlist_dir}\n")
        f.write(f"log_file_path: {log_file_path}\n")
        f.write("\n")
        f.write("Generating runlist...\n")
        f.write("\n")

    if runnumber > 0:
        # exicute the command to generate the runlist
        os.system(f'CreateDstList.pl --build {anabuild} -cdb {database} DST_CALOFITTING_run2auau --run {runnumber} > {log_file_path} 2>&1')
    else:
        # gen list of all runs
        os.system(f'CreateDstList.pl --build {anabuild} -cdb {database} DST_TRIGGERED_EVENT_run2auau --printruns > allruns.list')
        os.system(f'CreateDstList.pl --build {anabuild} -cdb {database} DST_TRIGGERED_EVENT_run2auau --list allruns.list > {log_file_path} 2>&1')
        os.system(f'rm allruns.list')

    # get a list of all files in the directory
    runlists = os.listdir(runlist_dir)

    # print the runlists to log file
    with open(log_file_path, "a") as f:
        f.write("\n")
        f.write("Runlists:\n")
        for runlist in runlists:
            f.write(runlist + "\n")
        f.write("\n")
    
    # split each runlist
    for runlist in runlists:
        file_path = os.path.join(runlist_dir, runlist)
        split_run_list(file_path, output_dir, n_file_per_list)
    
    # remove the runlist directory
    os.system(f'rm -rf {runlist_dir}')

    # cd back to the output directory
    os.chdir(output_dir)

    print (f"Runlist generation complete. Log file: {log_file_path}")
    return


if __name__ == "__main__":
    
    argparser = argparse.ArgumentParser()
    argparser.add_argument("-b", "--anabuild", help="anabuild", required=True)
    argparser.add_argument("-d", "--database", help="database", required=True)
    argparser.add_argument("-r", "--runnumber", help="run number", default=-1, type=int)
    args = argparser.parse_args()

    gen_run_list(args.anabuild, args.database, args.runnumber)

    sys.exit(0)


