#!/bin/ruby

=begin
Takashi Matsuda, 2014
MIT License	

Summing up values in csv files with a same parameter set.
=end

file=open(ARGV[0],'r')
puts "open :#{ARGV[0]}"
out_file=open(ARGV[1],'w')
puts "open : #{ARGV[1]}"
while line = file.gets
  vals=line.split(',')
  puts vals
  #for data containing the alpha row
  out_file.write("#{vals[0]},#{vals[1]},#{vals[2]},#{vals[3]},#{1-vals[4].to_f},#{vals[5]}")
end
file.close
out_file.close
puts "<<<<<FIN>>>>>"
