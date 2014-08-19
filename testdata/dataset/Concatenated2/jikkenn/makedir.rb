#! /usr/bin/ruby

files = Dir::entries(Dir::pwd)
for filepath in files
	if /id.fa$/ =~ filepath
		dirname = (filepath.split("ortholog", 2))[0]
		Dir::mkdir(dirname)
	end
end