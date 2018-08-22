source("./function.R")

raw_pH5_c234   <- "./curation/out/234/pH5_c234_159.fasta"
raw_pH5_c232   <- "./curation/out/232/pH5_c232_204.fasta"
raw_pH5Nx_c234 <- "./nx/out/pH5Nx_c234_240.fasta"
raw_pH5_c7     <- "./clade7/out/pH5Nx_c7_65.fasta" 

csv_eco_c234   <- "./state/out/eco_234.csv"
csv_eco_c232   <- "./state/out/eco_232.csv"


# stratification ------
# 234
sub_pH5_c234_159        <- taxaInfo( file = raw_pH5_c234, useTree = FALSE )
sub_pH5_c234_159[[ 8 ]] <- read.csv( csv_eco_c234, stringsAsFactors = FALSE )$states[ 
  match( sub_pH5_c234_159[[6]], read.csv( csv_eco_c234, stringsAsFactors = FALSE )$name ) ]
sub_pH5_c234_159[[ 8 ]] <- ifelse( startsWith( sub_pH5_c234_159[[ 8 ]], "D" ), "D", "W" )

ac_c234h5_D <- acSearch( faslist = sub_pH5_c234_159, keyword.dir = 8, keyword = "D")
ac_c234h5_W <- acSearch( faslist = sub_pH5_c234_159, keyword.dir = 8, keyword = "W")
subfastaSeq( AC = TRUE, AC_list = ac_c234h5_D, filedir = raw_pH5_c234 )
subfastaSeq( AC = TRUE, AC_list = ac_c234h5_W, filedir = raw_pH5_c234 )

# 232
sub_pH5_c232_204        <- taxaInfo( file = raw_pH5_c232, useTree = FALSE )
sub_pH5_c232_204[[ 8 ]] <- read.csv( csv_eco_c232, stringsAsFactors = FALSE )$states[ 
  match( sub_pH5_c232_204[[6]], read.csv( csv_eco_c232, stringsAsFactors = FALSE )$name ) ]
sub_pH5_c232_204[[ 8 ]] <- ifelse( startsWith( sub_pH5_c232_204[[ 8 ]], "D" ), "D", "W" )

ac_c232h5_D <- acSearch( faslist = sub_pH5_c232_204, keyword.dir = 8, keyword = "D", range = c( 2000, 2013 ), range.dir = 4 )
ac_c232h5_W <- acSearch( faslist = sub_pH5_c232_204, keyword.dir = 8, keyword = "W", range = c( 2000, 2013 ), range.dir = 4 )
subfastaSeq( AC = TRUE, AC_list = ac_c232h5_D, filedir = raw_pH5_c232 )
subfastaSeq( AC = TRUE, AC_list = ac_c232h5_W, filedir = raw_pH5_c232 )


# table ------
# 234

sub_pH5Nx_c234_240        <- taxaInfo( file = raw_pH5Nx_c234, useTree = FALSE )
sub_pH5Nx_c234_240[[ 8 ]] <- read.csv( csv_eco_c234, stringsAsFactors = FALSE )$states[ 
  match( sub_pH5Nx_c234_240[[6]], read.csv( csv_eco_c234, stringsAsFactors = FALSE )$name ) ]
sub_pH5Nx_c234_240[[ 8 ]] <- ifelse( startsWith( sub_pH5Nx_c234_240[[ 8 ]], "D" ), "D", "W" )
eco.234nx                 <- data.frame( id = sub_pH5Nx_c234_240[[ 6 ]], states = sub_pH5Nx_c234_240[[ 8 ]], stringsAsFactors = FALSE )
write.table( x = eco.234nx, file = "./beast/eco.234nx", sep = "\t", quote = FALSE, row.names = FALSE )


# 232
eco.232      <- data.frame( id = sub_pH5_c232_204[[ 6 ]], states = sub_pH5_c232_204[[ 8 ]], stringsAsFactors = FALSE )
eco.232_2014 <- eco.232[ which( floor( as.numeric( str_match( eco.232$id, "[0-9.]+$" ) ) ) < 2014 ), ]
write.table( x = eco.232_2014, file = "./beast/eco.232", sep = "\t", quote = FALSE, row.names = FALSE )


# 7
rmDup( raw_pH5_c7, year = c(2000, 2014) )
sub_pH5_c7_2014      <- taxaInfo( file = "./clade7/out/pH5nx_c7_57.fasta", useTree = FALSE )
sub_pH5_c7_2014[[8]] <- read.csv( "./state/out/eco_7.csv", stringsAsFactors = FALSE )$states[
  match( sub_pH5_c7_2014[[6]], read.csv( "./state/out//eco_7.csv", stringsAsFactors = FALSE )$name ) ]
sub_pH5_c7_2014[[8]] <- ifelse( startsWith( sub_pH5_c7_2014[[8]], "D" ), "D", "W" )
eco.7_2014           <- data.frame( id = sub_pH5_c7_2014[[6]], states = sub_pH5_c7_2014[[8]], stringsAsFactors = FALSE )
write.table( x = eco.7_2014, file = "./beast//eco.7.nx", sep = "\t", quote = FALSE, row.names = FALSE )


# tip data sampling ------

timeDice( fas.dir     = "./clade7/out/pH5nx_c7_57.fasta", 
          timetab.dir = "./others/out/time_ha_7326", 
          ecotab.dir  =  "./beast/eco.7.nx" )




