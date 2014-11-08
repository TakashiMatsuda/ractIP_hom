#!/usr/bin/python

"""
Takashi Matsuda, 2014
MIT License

Python script to sort csv file having N columns by 1st column value.
"""
import sys

filename = sys.argv[1]
fp = open(filename, 'r')
print "open : "+filename
"""
put a row into list of vectors with two values 
"""
res = []
flag = False
for line in fp:
	if flag==False:
		flag=True
		continue
	else:
		tmp = line.split(',')
		tmp_res=[]
		for vl in tmp:
			tmp_res.append(float(vl))
		res.append(tmp_res)

# sort that list by first value
sortedres = sorted(res, key=lambda x : x[int(sys.argv[2])])

op = open("sorted"+filename, 'w')
for line in sortedres:
	counter = 0
	for line_val in line:
		if counter==0:
			op.write(str(line_val))
		else:
			op.write(','+str(line_val))
			if counter == len(line)-1:
				op.write('\n')
		counter = counter + 1

op.close()
