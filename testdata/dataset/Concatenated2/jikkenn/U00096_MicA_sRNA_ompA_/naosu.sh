#! /bin/bash

for file in `ls *.fa`
do
	sed -i -e "s/T/U/g" ${file}
done
