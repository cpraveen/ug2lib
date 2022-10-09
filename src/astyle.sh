#!/bin/bash

files=`ls -1 *.cc *.h`

for f in $files
do
   echo $f
   astyle --options=./astyle.rc $f
done
