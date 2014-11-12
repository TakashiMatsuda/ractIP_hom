#! /usr/bin/ruby
# -*- coding: utf-8 -*-

=begin
Takashi Matsuda, 2014
MIT License
=end

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


# 回文構造
def innerbplist(structure)
  res = []
  nums=[]
  for i in 0..(structure.length-1)
    if structure[i] == '('
      nums << i
      next
    elsif structure[i] == ')'
      first = nums.pop
      res << [first, i]
      next
    end
  end
  return res
end


#非回文構造
def outerbplist(structure)
  fposlist = nums_instr('[', structure)
  bposlist = nums_instr(']', structure)
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
      if testbp[0] == ansbp[0]
        for ansbp2 in anslist
          if testbp[1] = ansbp2[1]
            tp = tp + 1
            break
          end
        end
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
  if /^ibp/ =~ filepath
  if /.fa$/ =~ filepath
#    puts "File Open: #{filepath}"
    test_file = open(filepath, 'r')
    #replacement of 'answer' to '_res'
    #tf_name = filepath.gsub('answer', '_res')
#    test_file = open("res/#{tf_name}", 'r')
    test_strc = secondstrc(test_file)
    ans_file = open("../MicA-ompAanswer.fa")
    ans_strc = secondstrc(ans_file)
    inbp = [innerbplist(test_strc), innerbplist(ans_strc)]# 2本の区別の情報は潰している。
    outbp = [outerbplist(test_strc), outerbplist(ans_strc)]
    # count the number of matching between bp[0](experiment) and bp[1](answer)
    inner_res=countmatch(inbp, ans_strc.length)
    outer_res=countmatch(outbp, ans_strc.length)
#tn, tp, fn and fp of innerbp
    output = open("#{filepath}-count-inner", 'w+')
    output.write("tn,#{inner_res[1]}\ntp,#{inner_res[0]}\nfn,#{inner_res[3]}\nfp,#{inner_res[2]}")
    output.close

    output = open("#{filepath}-count-total-sensitivity", 'w+')
    tp_tot = inner_res[0]+outer_res[0]
    fn_tot = inner_res[3]+outer_res[3]
#    puts "#{filepath}====tp: #{tp_tot}"
    sen_tot = Float(tp_tot) / Float(tp_tot + fn_tot)
    output.write("#{sen_tot}")
    output.close

    output = open("total-count-sen-ppv-fmeasure.csv", 'a+')
    # もしoutputが空だったら、先頭タグを書いてあげる。
    if output.read(1) == nil
      output.write("gamma_s,gamma_h,sensitivity,ppv,fmeasure\n")
    end
    # parameterの値、sen, ppv, fmeasureの順のcsvを作る
    fp_tot = inner_res[2] + outer_res[2]
    ppv_tot = Float(tp_tot) / Float(tp_tot + fp_tot)
    fmeasure=2*sen_tot*ppv_tot / (sen_tot + ppv_tot)

    # filepathをparseして、パラメータibp, hbp, alphaの値を得る
    parsedary = filepath.split('-')
    gamma_s = parsedary[1].to_f
    gamma_h = parsedary[3].to_f
# alphaを動かしている場合

    alpha=parsedary[5]
    alpha.slice!(-3..-1)
    output.write("#{gamma_s},#{gamma_h},#{alpha},#{sen_tot},#{ppv_tot},#{fmeasure}\n")

# alphaを動かしていない場合
=begin
    output.write("#{gamma_s},#{gamma_h},#{sen_tot},#{ppv_tot},#{fmeasure}\n")

=end

#    output.write("#{filepath},#{sen_tot},#{ppv_tot},#{fmeasure}\n")
    output.close

=begin
    output = open("plot.csv", 'a')
    tot = (test_strc.length * (test_strc.length - 1)) / 2
    output.write("#{inner_res[3]/tot},#{inner_res[1]/tot}\n")
    output.close
=end
#    puts "File close: #{filepath}"
#    tcount = tcount + tc_tmp
#    ncount = ncount + nc_tmp
#    puts "True : #{tc_tmp}"
#    puts "False: #{nc_tmp}"
  end
  end
end

puts "----------FIN-----------"
#puts "True: #{tcount}"
#puts "False: #{ncount}"
