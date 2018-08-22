source("./function.R")

# extract GsGD and rebuild another tree
ls.gsgd <- tagExtra( "./gsgd/in/pH5_7326_e0126.tre" )$id[ !is.na(tagExtra( "./gsgd/in/pH5_7326_e0126.tre" )$tag) ]
leafEx( "./data/out/trim_pH5_7326.fasta", ls.gsgd )


# clade annotation 
ann_6407 <- "./gsgd/out/pH5_7326_gsgd_e0702.tre"
pH5_seq  <- "./data/out/trim_pH5_7326.fasta"
pN1_seq  <- "./data/out/trim_pN1_4468.fasta"
p.table  <- "./data/out/pH5NA.csv"

ann_6407.tag     <- tagExtra( ann_6407 )
ann_6407.tag$tag <- gsub("ff0000", "c234", ann_6407.tag$tag)
ann_6407.tag$tag <- gsub("00ff00", "c232", ann_6407.tag$tag)
pTable           <- read.csv( p.table, stringsAsFactors = FALSE)


# clade 234
ac_all_c234 <- str_match( ann_6407.tag$id[ which( ann_6407.tag$tag == "c234" ) ], "[A-Z0-9]+" )[,1]
ac_h5_c234  <- pTable$ac.ha[ na.omit( match( ac_all_c234, pTable$ac.ha ) ) ]
ac_n1_c234  <- pTable$ac.na[ intersect( na.omit( match( ac_all_c234, pTable$ac.ha ) ), 
                                        which( pTable$sero == "H5N1") ) ]
subfastaSeq( AC = TRUE, filedir = pH5_seq, AC_list = ac_h5_c234 )
subfastaSeq( AC = TRUE, filedir = pN1_seq, AC_list = ac_n1_c234 )

# clade 232
ac_all_c232 <- str_match( ann_6407.tag$id[ which( ann_6407.tag$tag == "c232" ) ], "[A-Z0-9]+" )[,1]
ac_h5_c232  <- pTable$ac.ha[ na.omit( match( ac_all_c232, pTable$ac.ha ) ) ]
ac_n1_c232  <- pTable$ac.na[ intersect( na.omit( match( ac_all_c232, pTable$ac.ha ) ), 
                                        which( pTable$sero == "H5N1") ) ]
subfastaSeq( AC = TRUE, filedir = pH5_seq, AC_list = ac_h5_c232 )
subfastaSeq( AC = TRUE, filedir = pN1_seq, AC_list = ac_n1_c232 )


# rmdup 234
rmDup( fasfile = "./gsgd/out/pH5_c234_2429.fasta", rmdup = TRUE )
rmdup_plus( "./gsgd/out/pH5_c234_2429_cr.fasta" )

# rmdup 232
rmDup( fasfile = "./gsgd/out/pH5_c232_1296.fasta", rmdup = TRUE )
rmdup_plus( "./gsgd/out/pH5_c232_1296_cr.fasta" )

# for gsgd tree
rmdup_plus( "./gsgd/out/pH5_7326_gsgd.fasta" )