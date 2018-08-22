source("function.R")

require(seqinr)
require(stringr)
require(ape)
require(ggtree)
require(dplyr)


# clade annotation 
ann_6407 <- "./gsgd/out/pH5_7326_gsgd_e0702.tre"
pH5_seq  <- "./data/out/trim_pH5_7326.fasta"
pN1_seq  <- "./data/out/trim_pN1_4468.fasta"
p.table  <- "./data/out/pH5NA.csv"

ann_6407.tag     <- tagExtra( ann_6407 )
ann_6407.tag$tag <- gsub( "ff8000", "c7", ann_6407.tag$tag )
pTable           <- read.csv( p.table, stringsAsFactors = FALSE )

# clade 7
ac_all_c7 <- str_match( ann_6407.tag$id[ which( ann_6407.tag$tag == "c7" ) ], "^[A-Z0-9]+" )[,1]
ac_h5_c7  <- pTable$ac.ha[ na.omit( match( ac_all_c7, pTable$ac.ha ) ) ]
ac_n1_c7  <- pTable$ac.na[ intersect( na.omit( match( ac_all_c7, pTable$ac.ha ) ), 
                                      which( pTable$sero == "H5N1") ) ]

subfastaSeq( AC = TRUE, filedir = pH5_seq, AC_list = ac_h5_c7 )
subfastaSeq( AC = TRUE, filedir = pN1_seq, AC_list = ac_n1_c7 )

# rmdup 7
# rmDup( fasfile = "./clade7/pH5_c7_80.fasta", rmdup = TRUE )
# rmDup( fasfile = "./clade7/pN1_c7_63.fasta", rmdup = TRUE )


# curation --------

# 1 keep only China (+HK) seq.
# no obvious reintroduction case was observed

rmDup( fasfile = "./clade7/out//pH5_c7_80.fasta", sero = "H5N1", geo = c( "China", "Hong_Kong" ), rmdup = FALSE )
rmDup( fasfile = "./clade7/out//pN1_c7_63.fasta", sero = "H5N1", geo = c( "China", "Hong_Kong" ), rmdup = FALSE )

# 2 reassortment 
# no obvious reintroduction case was observed

# tre.pH5_7326_c7  <- read.nexus( "./clade7/out/pH5_c7.tre" )
# tre.pN1_4468_EA  <- read.nexus( "./curation/in/pN1_4468_EA.tre" )
# 
# cleantre.pH5_7326_c7           <-  tre.pH5_7326_c7
# cleantre.pH5_7326_c7$tip.label <-  sub( pattern = "^'[A-Z0-9]+_", replacement = "'", x = cleantre.pH5_7326_c7$tip.label )
# cleantre.pN1_4468_EA           <-  tre.pN1_4468_EA
# cleantre.pN1_4468_EA$tip.label <-  sub( "^'[A-Z0-9]+_", "'", cleantre.pN1_4468_EA$tip.label )
# 
# cleantxt_pH5       <- fortify( cleantre.pH5_7326_c7 )
# cleantxt_pH5       <- cleantxt_pH5[ which( cleantxt_pH5$isTip ), ]
# cleantxt_pH5$label <- gsub( "'", "", cleantxt_pH5$label )
# cleantxt_pH5       <- data.frame( cleantxt_pH5,  gene = "H5", stringsAsFactors = FALSE )
# 
# cleantxt_pN1       <- fortify( cleantre.pN1_4468_EA )
# cleantxt_pN1       <- cleantxt_pN1[ which( cleantxt_pN1$isTip ), ]
# cleantxt_pN1$label <- gsub( "'", "", cleantxt_pN1$label )
# cleantxt_pN1       <- data.frame( cleantxt_pN1,  gene = "N1", stringsAsFactors = FALSE )
# 
# cleantxt <- rbind( cleantxt_pH5, cleantxt_pN1 )
# cleantxt <- data.frame( cleantxt, geo  = str_match( cleantxt$label, "\\|([A-Za-z_]+)\\|" )[, 2],
#                         stringsAsFactors = FALSE )
# 
# cleantxt %>%
#   ggplot( aes(x = gene, y = y, group = label) ) + 
#   geom_point( size = 0.01 ) +
#   geom_line( alpha = 0.5 ) + 
#   theme_bw() + ylab("") + xlab("") + 
#   geom_hline( yintercept = 3150, color = "blue" ) +
#   geom_hline( yintercept = 2658, color = "blue" ) +
#   geom_hline( yintercept = 2100, color = "blue" ) +
#   scale_y_continuous( limits = c(0,300) )
# 
# tem.match        <- match( cleantxt_pH5$label, cleantxt_pN1$label )
# cleantxt_pNA_7   <- cleantxt_pN1[ tem.match[ !is.na( tem.match ) ], ]
# cleantxt_pHA_7   <- cleantxt_pH5[ which( cleantxt_pH5$label %in% cleantxt_pNA_7$label ), ]
# tem.list         <- intersect(  cleantxt_pHA_7$label, cleantxt_pNA_7$label )
# 
# ls.7   <- grep( "\\|China\\||\\|Hong_Kong\\|", tem.list, value = TRUE, ignore.case = TRUE )

# exactly as rmdup::h5n1 + cn_hk
# leafEx( filedir = "./clade7/out//pH5_c7_80.fasta", grep( paste0( gsub( "\\|", "\\\\|", ls.7 ), collapse = "|" ), fastaEx( "./clade7/out/pH5_c7_80.fasta" )$id, value = TRUE ) )
# leafEx( filedir = "./clade7/out/pN1_c7_63.fasta", grep( paste0( gsub( "\\|", "\\\\|", ls.7 ), collapse = "|" ), fastaEx( "./clade7/out/pN1_c7_63.fasta" )$id, value = TRUE ) )


# 3 rmDup
rmdup_plus( "./clade7/out/pH5_c7_80_52.fasta" ) #result = 50
rmdup_plus( "./clade7/out/pN1_c7_63_52.fasta" ) #result = 48       

# bigtree
rmdup_plus( "./clade7/out/pH5_c7_80.fasta")


# h5nx --------
# remove one strain of apparent outlier in the evolutionary trajectory 
rmDup( "./clade7/out//pH5_c7_80.fasta", geo = c( "China", "Hong_Kong" ), rmdup = FALSE)
rmdup_plus( "./clade7/out//pH5_c7_80_cr.fasta" )



