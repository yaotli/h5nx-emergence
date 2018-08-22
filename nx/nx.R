source("function.R")

require(seqinr)
require(stringr)

raw234h5 <- "./gsgd/out/pH5_c234_2429.fasta"
c2344tre <- "./nx/in/pH5_c234_2429_e0220.tre"

leafEx( raw234h5, tagExtra( c2344tre )$id[ !is.na( tagExtra( c2344tre )$tag ) ] )
rmDup( "./nx/out/pH5_c234_2429_1850.fasta", year = c(2000, 2014), geo = c("China", "Hong_Kong"), rmdup = TRUE )
