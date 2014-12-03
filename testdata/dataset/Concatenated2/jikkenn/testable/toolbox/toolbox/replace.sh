#!/bin/bash
sed -e s/$1/$2/g qsubjobedtest.sh > qsubjobedtest2.sh
rm qsubjobedtest.sh
mv qsubjobedtest2.sh qsubjobedtest.sh

