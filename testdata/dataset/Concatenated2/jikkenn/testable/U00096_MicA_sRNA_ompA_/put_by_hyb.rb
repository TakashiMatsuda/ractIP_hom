 #!/usr/bin/ruby

files = Dir::entries(Dir::pwd)
for filepath in files
  # find a answer file
  if /.fa$/ =~ filepath
    if /hyb-0.fa/ =~ filepath
      `mv #{filepath} ./hyb-0/`
    elsif /hyb-1.fa/ =~ filepath
      `mv #{filepath} ./hyb-1/`
  	elsif /hyb-0.1/ =~ filepath
  		`mv #{filepath} ./hyb-0.1/`
  	elsif /hyb-0.2/ =~ filepath
  		#statements
  		`mv #{filepath} ./hyb-0.2/`
  	elsif /hyb-0.3/ =~ filepath
  		#statements
  		`mv #{filepath} ./hyb-0.3/`
  	elsif /hyb-0.4/ =~ filepath 
  		#statements
  		`mv #{filepath} ./hyb-0.4/`
  	elsif /hyb-0.5/ =~ filepath
  		#statements
  		`mv #{filepath} ./hyb-0.5/`
  	elsif /hyb-0.6/ =~ filepath
  		#statements
  		`mv #{filepath} ./hyb-0.6/`
  	elsif /hyb-0.7/ =~ filepath
  		#statements
  		`mv #{filepath} ./hyb-0.7/`
  	elsif /hyb-0.8/ =~ filepath
  		#statements
  		`mv #{filepath} ./hyb-0.8/`
  	elsif /hyb-0.9/ =~ filepath
  		#statements
  		`mv #{filepath} ./hyb-0.9/`
  	elsif /hyb-1.0/ =~ filepath
  		#statements
  		`mv #{filepath} ./hyb-1.0/`

    end
#    test_file = open(filepath, 'r')
	# mv file to dir by hyb-w
  end
end

puts "fin"
