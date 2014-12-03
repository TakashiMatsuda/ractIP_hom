#!/bin/ruby

=begin
Takashi Matsuda, 2014
MIT License	

Summing up values in csv files with a same parameter set.
=end

micA_ompA_file=open("total-count-sen-ppv-fmeasure-MicA-ompA.csv", 'r')
oxyS_fhlA_file=open("total-count-sen-ppv-fmeasure-OxyS-fhlA.csv", 'r')
ryhB_sodB_file=open("total-count-sen-ppv-fmeasure-RyhB-sodB.csv", 'r')
averaged_file=\
	open("total-count-sen-ppv-fmeasure-averaged.csv", 'w')
mA = micA_ompA_file.gets
oS=oxyS_fhlA_file.gets
rB=ryhB_sodB_file.gets
averaged_file.write(mA)
while mA = micA_ompA_file.gets
  oS=oxyS_fhlA_file.gets
  rB=ryhB_sodB_file.gets
  mA_vals=mA.split(',')
  oS_vals=oS.split(',')
  rB_vals=rB.split(',')
  if not(mA_vals[0]==oS_vals[0] and oS_vals[0]==rB_vals[0] \
         and mA_vals[1]==oS_vals[1] and oS_vals[1]==rB_vals[1])
    print "ALERT:: 3 files do not contain same parameters"
  end
  averaged_file.write("#{mA_vals[0]},#{mA_vals[1]},#{(mA_vals[2]+oS_vals[2]+rB_vals[2])/3.0},#{(mA_vals[3]+oS_vals[3]+rB_vals[3])/3.0},#{(mA_vals[4]+oS_vals[4]+rB_vals[4])/3.0}\n")
end
