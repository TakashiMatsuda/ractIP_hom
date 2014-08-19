#! /usr/bin/ruby
require 'open4'

files = Dir::entries(Dir::pwd)
for filepath in files
	if /fa$/ =~ filepath
		puts filepath
		#{}`/Users/takashi/cbrc/ractip_hom/testdata/dataset/sim2.py #{filepath}`
		namelist = []
		tmp_split = []
		tmp_split = filepath.split(/-/)
		namelist << tmp_split[0]
		namelist << ((tmp_split[1]).split(/answer/))[0]
		paramlist = ["alpha", "beta", "mix-weight", "hyb-mix-weight"]
		for param in paramlist
			for i in 0..10
				result = `/Users/takashi/cbrc/ractip_hom/src/ractip_hom --#{param}=#{i/10} tmp/#{namelist[0]}.fa tmp/#{namelist[0]}_simulation.fa tmp/#{namelist[1]}.fa tmp/#{namelist[1]}_simulation.fa`
				puts result
				system("mkdir res")
				res_file = open("res/#{namelist[0]}-#{namelist[1]}_res_alpha_#{i/10}.fa", 'w')
				res_file.write(result)
				res_file.close()
			end
		end
	end
end
puts "fin"