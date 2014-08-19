#!/bin/bash

files="./*"
mkdir useful
for filepath in ${files}
do
	echo ${filepath}
	tr '&' '\n' < ${filepath}
done
