 #!/usr/bin/ruby

files = Dir::entries(Dir::pwd)
for filepath in files
  # find a answer file
  if /.fa$/ =~ filepath
    if /hyb-0.fa/ =~ filepath
      `mkdir hyb-0`
      `mv #{filepath} ./hyb-0/`
    elsif /hyb-1.fa/ =~ filepath
      `mkdir hyb-1`
      `mv #{filepath} ./hyb-1/`
  	elsif /hyb-0.1/ =~ filepath
      `mkdir hyb-0.1`
  		`mv #{filepath} ./hyb-0.1/`
  	elsif /hyb-0.2/ =~ filepath
  		#statements
      `mkdir hyb-0.2`
  		`mv #{filepath} ./hyb-0.2/`
  	elsif /hyb-0.3/ =~ filepath
  		#statements
      `mkdir hyb-0.3`
  		`mv #{filepath} ./hyb-0.3/`
  	elsif /hyb-0.4/ =~ filepath 
  		#statements
      `mkdir hyb-0.4`
  		`mv #{filepath} ./hyb-0.4/`
  	elsif /hyb-0.5/ =~ filepath
  		#statements
      `mkdir hyb-0.5`
  		`mv #{filepath} ./hyb-0.5/`
  	elsif /hyb-0.6/ =~ filepath
  		#statements
      `mkdir hyb-0.6`
  		`mv #{filepath} ./hyb-0.6/`
  	elsif /hyb-0.7/ =~ filepath
  		#statements
      `mkdir hyb-0.7`
  		`mv #{filepath} ./hyb-0.7/`
  	elsif /hyb-0.8/ =~ filepath
  		#statements
      `mkdir hyb-0.8`
  		`mv #{filepath} ./hyb-0.8/`
  	elsif /hyb-0.9/ =~ filepath
  		#statements
      `mkdir hyb-0.9`
  		`mv #{filepath} ./hyb-0.9/`
  	elsif /hyb-1.0/ =~ filepath
  		#statements
      `mkdir hyb-1.0`
  		`mv #{filepath} ./hyb-1.0/`

    end
#    test_file = open(filepath, 'r')
	# mv file to dir by hyb-w
  end
end

puts "fin"
