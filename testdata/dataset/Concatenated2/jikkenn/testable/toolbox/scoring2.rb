#! /usr/bin/ruby

# f: filepointer
# pick a RNA secondstructure section from fasta file
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


def nums_instr(ch, str)
  nums = []
  for i in 0..(str.length-1)
    if str[i] == ch
      nums << i
    end
  end
  return nums
end


def innerbplist(structure)
  fposlist = nums_instr('(', structure)
  bposlist = nums_instr(')', structure)
  res = []
  for i in 0..(fposlist.length-1)
    res << [fposlist[i], bposlist[i]]
  end
  return res
end


def countmatch(bplist_list, l)
  testlist = bplist_list[0]
  anslist = bplist_list[1]
  tp = 0
  fp = 0
  tn = 0
  fn = 0
  for testbp in testlist
    for ansbp in anslist
      if testbp[0] < ansbp[0]
        next
      elsif testbp[0] == ansbp[0]
        for ansbp2 in anslist
          if testbp[1] < ansbp2[1]
            next
          elsif testbp[1] = ansbp2[1]
            tp = tp + 1
            break
          elsif testbp[1] > ansbp2[1]
#            fp = fp + 1
            break
          end      
        end
        break
      elsif testbp[0] > ansbp[0]
#        fp = fp + 1
        break
      end
    end
  end
  fp = testlist.length - tp
  fn = anslist.length - tp
  tn = (l*(l-1))/2 - tp - fp - fn
  return [tp, tn, fp, fn]
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
    test_strc = secondstrc(test_file)
    ans_file = open("../../OxyS-fhlAanswer.fa", 'r')
    ans_strc = secondstrc(ans_file)

    puts ans_strc
    puts test_strc
    inbp = [innerbplist(test_strc), innerbplist(ans_strc)]# 2本の区別の情報は潰している。
    outbp = [outerbplist(test_strc), outerbplist(ans_strc)]
    # count the number of matching between bp[0](experiment) and bp[1](answer)
    inner_res=countmatch(inbp)
    outer_res=countmatch(outbp)

=begin
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
=end
    output = open("#{filepath}-count-inner", 'w+')
    output.write("tn,#{inner_res[1]}\ntp,#{inner_res[1]}\nfn,#{inner_res[3]}\nfp,#{inner_res[2]}")
    output.close
    output = open("plot.csv", 'a')
    tot = (test_strc.length * (test_strc.length - 1)) / 2
    output.write("#{inner_res[3]/tot},#{inner_res[1]/tot}\n")
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
