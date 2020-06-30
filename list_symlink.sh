#!/bin/bash

export LANG=C

cd `dirname $0`
SCRIPT_DIR=`pwd`

cd $SCRIPT_DIR

rm -f symlink_list.txt

while read line
do
  f=`cygpath $line`
  f=`echo $f | sed -e "s/[\r\n]\+//g"`
  f_=`readlink $f`

  echo $f $f_ >> symlink_list.txt
  rm -f $f

done < junction_list.txt

cp make_symlink.cmd make_symlink.sh symlink_list.txt /

# Remove garbage files
rm -rf /home/*

date
