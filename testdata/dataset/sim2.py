#! /usr/bin/python
#! -*- coding:utf-8 -*-

import re 
import random
import sys
import commands

N = 20
P = 0.4
T = 10000


## seq_list[0],set_list[1]に配列が格納された
## それぞれに対して変異データを確率的に起こす。サンプリングを行う。
def variationmaker(n, p, t, seq, RD):
	## nは要求するアラインメントに含まれる配列の本数
	## pは単塩基変異が単位時間あたりに起こる確率
	## tは経過する時間
	if RD == "RNA":
		pyrimidine='U'
	else:
		pyrimidine='T'

	alignment = []
	##乱数を発生
	NN = 1000000
	random.seed()
	for i in range(n):
		print str(i)+'/'+str(n)
		str_changed=""
		for j in range(t):
			for count, s in enumerate(seq):
				r = random.random()
				## 以下の書き方を変更して、意味のあるものにする。
				if (p/3) > r:
					if s=='A':
						str_changed = seq[:count] + 'C' + seq[count+1 :]
					elif s == 'C' and s == 'G' and s == pyrimidine:
						str_changed = seq[:count] + 'A' + seq[count+1 :]
				elif (p/3) <= r and r < (p*2/3):
					if s == 'C':
						str_changed = seq[:count] + 'G' + seq[count+1 :]
					elif s == 'A':
						str_changed = seq[:count] + 'G' + seq[count+1 :]
					elif s == 'G' and s == pyrimidine:
						str_changed = seq[:count] + 'C' + seq[count+1 :]
				elif (p*2/3) <= r and r < p :
					if s == 'A':
						str_changed = seq[:count] + pyrimidine + seq[count+1 :]
					elif s == 'C':
						str_changed = seq[:count] + pyrimidine + seq[count+1 :]
					elif s == 'G':
						str_changed = seq[:count] + pyrimidine + seq[count+1 :]
					elif s == pyrimidine:
						str_changed = seq[:count] + 'G' + seq[count+1 :]
		alignment.append(str_changed)
		print str_changed
	return alignment


sys.argv[1]


file_name = sys.argv[1]
seq_name = []
p = re.compile('-')
tmp_list = p.split(file_name)
seq_name.append(tmp_list[0])
p = re.compile('answer')
seq_name.append(p.split(tmp_list[1])[0])

seq_file=open(file_name,'r')
nametag_pattern = '^>'
seq_pattern = '^A|^U|^G|^C|^N'
nametag = []
seq_list = []
for line in seq_file:
	## matching
	match_tag = re.search(nametag_pattern, line)
	if match_tag:
		nametag.append(line)
	else:
		match_seq = re.search(seq_pattern, line)
		if match_seq:
			seq_list.append(line)

seq_file.close()
alignments=[]

alignments.append(variationmaker(N, P, T, seq_list[0], "RNA"))
alignments.append(variationmaker(N, P, T, seq_list[1], "RNA"))

## write operation
L = len(alignments)
i=0
for align in alignments:
	commands.getoutput("mkdir tmp")
	output_file = open("tmp/"+seq_name[i]+".fa", 'w')
	output_file.write(nametag[i])
	output_file.write(seq_list[i])
	output_file.close()

	output_file = open("tmp/"+seq_name[i]+"_simulation.fa", 'w')
	output_file.write(nametag[i])
	output_file.write(seq_list[i])
	for count, s in enumerate(align) :
		output_file.write(nametag[i][:len(seq_name[i])-1]+"_sim_"+str(count)+'\n')
		output_file.write(s)
	output_file.close()
	i += 1
