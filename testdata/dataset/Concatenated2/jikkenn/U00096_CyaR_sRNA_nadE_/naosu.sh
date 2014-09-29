#! /bin/bash

for file in `ls *.fa`
do
	sed -i -e "s/-/./g" ${file}
done
