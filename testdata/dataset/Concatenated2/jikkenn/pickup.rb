 #! /usr/bin/ruby

 # アラインメントfastaから順番に一つづつ各配列を抜き出す。

filename = ARGV[0]
file_point = open(ARGV[0], 'r')
content = file_point.read()

seqlist = content.split('>')

