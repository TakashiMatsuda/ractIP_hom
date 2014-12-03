#!/bin/sh

mkdir tritest_hyb_t_u
for i in 0.0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.55 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.4 2 3 4 5 8 10
do
    echo "alpha == ${i}"
    qsub  qsubjobedtest.sh ${i}
done
