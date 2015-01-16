#!/bin/bash

if [ $# -ne 1 ] 
then
	echo "usage: ./ezfslbin.sh bin_dir_to_create_or_use"
	exit
else
	fslbin=$1
fi

echo "this will make the directory $fslbin will many symbolic links inside"
echo "sleeping 2 seconds before executing..."

sleep 2

echo "mkdir $fslbin"
set -o verbose
mkdir -p $fslbin


known_fsl_cmd=eddy_correct
echo "assuming fsl cmd $known_fsl_cmd will be found"
known_fsl_exe=fsl5.0-$known_fsl_cmd
echo "assuming $known_fsl_exe exists..."
path_to_known=$(which $known_fsl_exe)
echo "found $path_to_known"
base_pattern=${path_to_known%$known_fsl_cmd}
echo "using base pattern $base_pattern"
for f in $(ls $base_pattern*);
do
	target=$(pwd)/$fslbin/${f#$base_pattern}
	echo "creating link from $target to $f"
	ln -s $f $target
done


echo "do export PATH=\$(pwd)/$fslbin:\$PATH to use"
