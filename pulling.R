library(rgl)
data <- read.table("/Users/takashi/cbrc/ractip_hom/testdata/dataset/toysample/bp-matrix-wh-0.5/out_hp_2.csv", skip=2, sep=",")
mat <- matrix(1:20, ncol=4, nrow=5)
x <- c(1:ncol(data))
y <- c(1:nrow(data))
x
y

plot3d(x, y, mat)
scatterplot3d(x,y,data)
mode(mat)
mode(data)
la