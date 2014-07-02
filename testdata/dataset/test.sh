#!/bin/sh
files = "./*"
for filepath in ${files}
do
	./sim2.py ${filepath}
	../../../src/ractip_hom 