#! /usr/bin/ruby

# makes a side effect.副作用あり。
# '&' と'>'で分割し、それよりを前を抽出しfilename[0], [1]にに書き込むことで、
# RNA joint structure fasta file を２つのRNAに分割する。
def split_fasta_write_recursive(str, filename)

	# ２次構造の部分が流れてきたらfilename[3], filename[4]宛に格納する。
	tag = str.split("\n", 2)[0]
	if /interaction/ =~ tag
		splited_and = str.split('&')

		# 一本目の配列の構造の書き込み
		answer= open(filename[2], 'w')
		answer.write(">#{tag}\n")
		answer.write(splited_and[0])

		# 二本目の配列の構造の書き込み
		answer.write(">#{tag}\n")# たぶんこの行は入れる必要がない。
		answer.write(splited_and[1])
		answer.close()

		return 0
	else

		out = [open(filename[0], 'w'), open(filename[1], 'w')]
		if str.include?("&")
			splited_and = str.split('&', 2)
     		# 一行目はタグ
     		former = splited_and[0].split("\n", 2)
    		# open the former file
	    	out[0].write(">#{former[0]}\n") # write tag
		    out[0].write(former[1]) # write the former seq
		    out[0].close()

	    	# pick the latter up
     		if splited_and[1].include?(">")
     			latter_seq = str.split('>', 2)
     			out[1].write(">#{former[0]}\n")
     			out[1].write(latter_seq[0])
     			out[1].close()

     			split_fasta_write_recursive(latter[1], filename)
     		else
     			raise "unsuitable format"
     		end
     	end
     end
 end


files = Dir::entries(Dir::pwd)
for filepath in files
	if /id.fa$/ =~ filepath

		puts "File Open: #{filepath}"

		fpoint = open(filepath, 'r')
		content = fpoint.read()
		fpoint.close()

		dirname = (filepath.split("ortholog", 2))[0]
		Dir::mkdir(dirname)

		outputname = Array.new
		outputname.push("#{dirname}1.fa")
		outputname.push("#{dirname}2.fa")
		outputname.push("#{dirname}answer")

		split_fasta_write_recursive(content, outputname)

	end
end