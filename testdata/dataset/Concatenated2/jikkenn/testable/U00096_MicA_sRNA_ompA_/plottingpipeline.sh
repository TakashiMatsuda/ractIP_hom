#!/bin/bash
rm ./*-count
rm ./*.csv
ruby put_by_hyb.rb
ruby ../../scoring2.rb
python ../../sortcsv2.py plot.csv
python ../../plotting.py sortedplot.csv
