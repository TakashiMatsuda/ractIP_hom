#! /usr/bin/ruby

def secondstrc(f)
	c = 0
	res = ""
	while l = f.gets
		if (c % 3 == 2)
			res << l
		end
		c = c + 1
	end
	return res
end

files = Dir::entries(Dir::pwd)
tcount = 0
ncount = 0
tot = 0
for filepath in files
	if /answer.fa$/ =~ filepath
		puts "File Open: #{filepath}"
		ans_file = open(filepath, 'r')
		#replacement of 'answer' to '_res'
		tf_name = filepath.gsub('answer', '_res')
		test_file = open("res/#{tf_name}", 'r')
		ans_strc = secondstrc(ans_file)
		test_strc = secondstrc(test_file)
		puts ans_strc
		puts test_strc
		len= ans_strc.length
		tot = tot + len
		tc_tmp = 0
		nc_tmp = 0
		for i in 0..(len-1)
			if ans_strc[i] == test_strc[i] then
				tc_tmp = tc_tmp + 1
			else
				nc_tmp = nc_tmp + 1
			end
		end
		tcount = tcount + tc_tmp
		ncount = ncount + nc_tmp
		puts "True : #{tc_tmp}"
		puts "False: #{nc_tmp}"
	end
end

puts "----------RESULT-----------e"
puts "True: #{tcount}"
puts "False: #{ncount}"