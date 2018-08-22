source("function.R")

require(stringr)
require(seqinr)

raw_pH5_c232   <- "./curation/out/232/pH5_c232_204.fasta"
raw_pN1_c232   <- "./curation/out/232/pN1_c232_189.fasta"

rmDup( fasfile = raw_pH5_c232, year = c( 2000, 2014 ), rmdup = FALSE )
rmDup( fasfile = raw_pN1_c232, year = c( 2000, 2014 ), rmdup = FALSE )

# tip date sampling with timeDice, for example

timeDice( fas.dir     = "./curation/out/234/pH5_c234_159.fasta", 
          timetab.dir = "./others/out/time_ha_7326", ecotable = FALSE)



