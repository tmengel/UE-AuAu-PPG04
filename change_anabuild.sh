#!/bin/bash
currdir=$(pwd)
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <anabuild>"
    exit 1
fi


cleanInstall $INSTALLDIR $1
echo "New install directory: $INSTALLDIR"
echo "New ana build: $BUILDVER"
srcLocal 

cd $currdir/src

#get list of subdirs
subdirs=(eventselector underlyingevent ppg04base)
for dir in "${subdirs[@]}"; do
    echo "Building package: $dir"
    cd $currdir/src/$dir
    buildModule -n
    echo "done building package: $dir"
done
cd $currdir
echo "Building all packages done"
srcLocal