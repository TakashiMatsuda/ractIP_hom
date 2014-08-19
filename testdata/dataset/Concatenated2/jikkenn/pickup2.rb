#! /usr/bin/ruby

#filename = ARGV[0]

def pickup(filename, dirname)
	file = open("#{dirname}/#{filename}", 'r')

	if /1.fa/ =~ filename
		n = 1
	elsif /2.fa/ =~ filename
		n = 2
	else
		return 0
	end

	output = open("/Users/takashi/play.txt", 'r')
	while line = file.gets
		if /^>/ =~ line
			output.close()
			output = open("#{dirname}/run/#{line[1..-2]}_#{n}.fa", 'w')
		end
		output.write(line)
	end
	output.close()
	file.close()
end

dirs = Dir::entries(Dir::pwd)
for dirname in dirs
	if File::ftype(dirname) == "directory"
		if /^[[:upper:]]/ =~ dirname
			Dir::mkdir("#{dirname}/run")
			pickup("#{dirname}1.fa", dirname)
			pickup("#{dirname}2.fa", dirname)
		end
	end
end