#!/bin/bash

cd `dirname $0`

while read line
do
  set ${line}
  ln -s ${2} ${1}
done < symlink_list.txt

