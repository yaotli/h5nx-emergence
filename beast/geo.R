source("function.R")

# input --------

rawfas.234 <- "./gsgd/out/pH5_c234_2429.fasta"
rawfas.232 <- "./gsgd/out/pH5_c232_1296.fasta"
rawfas.7   <- "./clade7/out/pH5_c7_80.fasta"


# light cleaning --------

rmDup( rawfas.234, year = c( 1996, 3000), rmdup = TRUE )
rmDup( rawfas.232, year = c( 1996, 3000), rmdup = TRUE )
rmDup( rawfas.7, year = c( 1996, 3000), rmdup = TRUE )

rmdup_plus( "./gsgd/out/pH5_c234_2429_cr.fasta" )
rmdup_plus( "./gsgd/out/pH5_c232_1296_cr.fasta" )
rmdup_plus( "./clade7/out/pH5_c7_80_cr.fasta" )


# geo dist. --------

faslist_234 <- taxaInfo( "./gsgd/out/pH5_c234_2429_cr_rmdP.fasta" )
faslist_232 <- taxaInfo( "./gsgd/out/pH5_c232_1296_cr_rmdP.fasta" )
faslist_7   <- taxaInfo( "./clade7/out/pH5_c7_80_cr_rmdP.fasta" )

faslist_234[[8]] <- geoID( faslist_234[[2]] )
faslist_232[[8]] <- geoID( faslist_232[[2]] )
faslist_7[[8]]   <- geoID( faslist_7[[2]] )


# remove 1 from 232 (the only north americas isolate)


# subsampling with geo. as distinct property --------

# clade 7
cladeSampling( trefile   = "./clade7/out/pH5_c7_80_cr_rmdP.tre", 
               fasfile   = "./clade7/out/pH5_c7_80_cr_rmdP.fasta", 
               suppList  = TRUE, 
               listinput = faslist_7,
               grid      = 0.5, 
               list.x    = c(6,4,8), 
               saveFasta = TRUE )
# clade 234
cladeSampling( trefile   = "./gsgd/out/pH5_c234_2429_cr_rmdP.tre", 
               fasfile   = "./gsgd/out/pH5_c234_2429_cr_rmdP.fasta", 
               suppList  = TRUE, 
               listinput = faslist_234,
               grid      = 0.5, 
               list.x    = c(6,4,8), 
               saveFasta = TRUE )
# clade 232
cladeSampling( trefile   = "./gsgd/out/pH5_c232_1296_cr_rmdPe.tre", 
               fasfile   = "./gsgd/out/pH5_c232_1296_cr_rmdP.fasta", 
               suppList  = TRUE, 
               listinput = faslist_232,
               grid      = 0.5, 
               list.x    = c(6,4,8), 
               saveFasta = TRUE )


# prepare BEAST --------

# clade 7
geo.big.7 <- data.frame( id     = fastaEx( "./clade7/out/pH5_c7_80_cr_rmdP_s.fasta" )$id, 
                         states = faslist_7[[8]][ match( fastaEx( "./clade7/out/pH5_c7_80_cr_rmdP_s.fasta" )$id, faslist_7[[6]] ) ], 
                         stringsAsFactors = FALSE )

geo.big.7$states <- gsub( "^g", "", geo.big.7$states )
write.table( x = geo.big.7, file = "./beast/geo.7", sep = "\t", quote = FALSE, row.names = FALSE )

timeDice( "./clade7/out/pH5_c7_80_cr_rmdP_s.fasta", "./beast/geo.7", "./others/out/time_ha_7326" )


# clade 234
geo.big.234 <- data.frame( id     = fastaEx( "./gsgd/out/pH5_c234_2429_cr_rmdP_s.fasta" )$id, 
                           states = faslist_234[[8]][ match( fastaEx( "./gsgd/out/pH5_c234_2429_cr_rmdP_s.fasta" )$id, faslist_234[[6]] ) ], 
                           stringsAsFactors = FALSE )

write.table( x = geo.big.234, file = "./beast/geo.234", sep = "\t", quote = FALSE, row.names = FALSE )
timeDice( "./gsgd/out/pH5_c234_2429_cr_rmdP_s.fasta", "./beast/geo.234", "./others/out/time_ha_7326" )


# clade 232
geo.big.232 <- data.frame( id     = fastaEx( "./gsgd/out/pH5_c232_1296_cr_rmdP_s.fasta" )$id, 
                           states = faslist_232[[8]][ match( fastaEx( "./gsgd/out/pH5_c232_1296_cr_rmdP_s.fasta" )$id, faslist_232[[6]] ) ], 
                           stringsAsFactors = FALSE )

write.table( x = geo.big.232, file = "./beast/geo.232", sep = "\t", quote = FALSE, row.names = FALSE )
timeDice( "./gsgd/out/pH5_c232_1296_cr_rmdP_s.fasta", "./beast/geo.232", "./others/out/time_ha_7326" )





