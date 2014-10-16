#! /usr/bin/ruby

# f: filepointer
# pick a RNA secondstructure
def secondstrc(f)
  res = ""
  flag = 0
  while l = f.gets
    if flag == 2 then
      res << l
      flag = 0
    elsif flag == 1 then
      flag = 2
    elsif /^>/ =~ l then
      flag = 1
    end
  end
  return res
end



# Browse all files in the current directory
files = Dir::entries(Dir::pwd)
#tcount = 0
#ncount = 0
#tot = 0
for filepath in files
  # find a answer file
  if /.fa$/ =~ filepath
    puts "File Open: #{filepath}"
    test_file = open(filepath, 'r')
    #replacement of 'answer' to '_res'
    #tf_name = filepath.gsub('answer', '_res')
#    test_file = open("res/#{tf_name}", 'r')
    ans_file = open("../../OxyS-fhlAanswer.fa", 'r')
    ans_strc = secondstrc(ans_file)
    test_strc = secondstrc(test_file)
    puts ans_strc
    puts test_strc
    len= ans_strc.length
#    tot = tot + len
 #   tc_tmp = 0
#    nc_tmp = 0
    tn = 0
    tp = 0
    fn = 0
    fp = 0
    for i in 0..(len-1)
      if ans_strc[i] == test_strc[i] then
        # to add TP and TN
        if test_strc[i] == '.' then
          tn = tn + 1
        else
          tp = tp + 1
        end
#        tc_tmp = tc_tmp + 1
      else
        if test_strc[i] == '.' then
          fn = fn + 1
        else
          fp = fp + 1
#         nc_tmp = nc_tmp + 1
        end
      end
    end
    output = open("#{filepath}-count", 'w+')
    output.write("tn,#{tn}\ntp,#{tp}\nfn,#{fn}\nfp,#{fp}")
    output.close
    output = open("plot.csv", 'a')
    tot = Float(len)
    output.write("#{fp/tot},#{tp/tot}\n")
    output.close
#    tcount = tcount + tc_tmp
#    ncount = ncount + nc_tmp
#    puts "True : #{tc_tmp}"
#    puts "False: #{nc_tmp}"
  end
end

puts "----------FIN-----------"
#puts "True: #{tcount}"
#puts "False: #{ncount}"
