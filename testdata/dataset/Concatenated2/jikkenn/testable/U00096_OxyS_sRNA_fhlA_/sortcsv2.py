#!/usr/bin/python

"""
Takashi Matsuda, 2014
MIT License

Python script to sort csv file having 2 column by 1st column value.
"""
import sys

filename = sys.argv[1]
fp = open(filename, 'r')
print "open : "+filename
"""
put a row into list of vectors with two values 
"""
res = []
for line in fp:
	tmp = line.split(',')
	tmp_res=[float(tmp[0]), float(tmp[1])]
	res.append(tmp_res)

# sort that list by first value
sortedres = sorted(res, key=lambda x : x[0])

op = open("sorted"+filename, 'w')
for line in sortedres:
	op.write(str(line[0])+','+str(line[1])+'\n')

op.close()
