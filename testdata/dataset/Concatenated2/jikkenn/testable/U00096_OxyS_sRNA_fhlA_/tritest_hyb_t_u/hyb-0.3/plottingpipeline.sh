#!/bin/bash
ruby ../../scoring2.rb
python ../../sortcsv2.py plot.csv
python ../../plotting.py sortedplot.csv
