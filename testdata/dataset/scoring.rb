#! /usr/bin/ruby
# 一致する塩基対の数を表示するスクリプト.

def secondstrc(f)
	c = 0
	res = []
	while l = f.gets
		if (c % 3 == 2)
			res << l
		end
		c = c + 1
	end
	return res
end

#strc: RNA second structure (string)
def innerbp(strc)
	bp_stack = []
	bp_list = []

	len = strc.length
	for i in 0..(len-1)
		for j in 0..(strc[i].size-1)
			if strc[i].eql?("(")
				bp_stack << i
			elsif strc[i].eql?(")")
				bp_list << [bp_stack.delete_at(-1) , i]
			end
		end
	end

	return bp_list
end


def hybridbp(strc)
	nandemo = []
	bp_list = []

	for i in 0..(strc[0].size - 1)
		if strc[0][i] == "["
			nandemo << i
		end
	end
	for i in 0..(strc[1].size - 1)
		if strc[1][i] == "]"
			bp_list << [nandemo.delete_at(0), i]
		end
	end

	print bp_list

	return bp_list
end

files = Dir::entries(Dir::pwd)
tcount = 0
ncount = 0
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

		innerbp_ans_l = Object.new
		innerbp_test_l = Object.new
		innerbp_ans_l = innerbp(ans_strc)
		innerbp_test_l = innerbp(test_strc)

		hybrid_ans_l = Object.new
		hybrid_test_l = Object.new
		hybrid_ans_l = hybridbp(ans_strc)
		hybrid_test_l = hybridbp(test_strc)

		tc_tmp = 0
		nc_tmp = 0

		
		for i in 0..(innerbp_ans_l.size - 1)

			tc_tmp = tc_tmp + innerbp_test_l.count(innerbp_ans_l[i])
			tc_tmp = tc_tmp + hybrid_test_l.count(hybrid_ans_l[i])
		end
		tcount = tcount + tc_tmp
		#ncount = ncount + nc_tmp
		puts "True : #{tc_tmp}"
		#puts "False: #{nc_tmp}"
	end
end

puts "----------RESULT-----------e"
puts "True: #{tcount}"
puts "False: #{ncount}"
