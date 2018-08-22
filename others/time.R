source("function.R")
require(stringr)

# raw data

raw.id_G_h  <- fastaEx("data/in/pH5_G_2136_20171124.fasta")$id
t.raw_G_h   <- gsub( pattern = "_\\(Day_unknown\\)_", "-99_", raw.id_G_h)
t.raw_G_h   <- gsub( pattern = "_\\(Month_and_day_unknown\\)_", "-99-99_", t.raw_G_h)
t.raw_G_h   <- str_match( t.raw_G_h, "_([0-9]{4}-[0-9]{2}-[0-9]{2})_E")[,2]

raw.id_G_n  <- fastaEx("data/in/pNA_G_2135_20171124.fasta")$id
t.raw_G_n   <- gsub( pattern = "_\\(Day_unknown\\)_", "-99_", raw.id_G_n)
t.raw_G_n   <- gsub( pattern = "_\\(Month_and_day_unknown\\)_", "-99-99_", t.raw_G_n)
t.raw_G_n   <- str_match( t.raw_G_n, "_([0-9]{4}-[0-9]{2}-[0-9]{2})_E")[,2]

raw.id_N_h  <- fastaEx("data/in/pH5_N_5753_20171124.fasta")$id
t.raw_N_h   <- gsub( pattern = "--$", "-99-99", raw.id_N_h)
t.raw_N_h   <- gsub( pattern = "-$", "-99", t.raw_N_h)
t.raw_N_h   <- str_match( t.raw_N_h, "_([0-9]{4}-[0-9]{2}-[0-9]{2})$")[,2]

raw.id_N_n  <- fastaEx("data/in/pNA_N_5721_20171124.fasta")$id
t.raw_N_n   <- gsub( pattern = "--$", "-99-99", raw.id_N_n)
t.raw_N_n   <- gsub( pattern = "-$", "-99", t.raw_N_n)
t.raw_N_n   <- str_match( t.raw_N_n, "_([0-9]{4}-[0-9]{2}-[0-9]{2})$")[,2]

t.raw_h  <- c( t.raw_G_h, t.raw_N_h )
t.raw_n  <- c( t.raw_G_n, t.raw_N_n )

ac.raw_h <- c( gsub( "EPI_ISL_", "EPI", str_match( raw.id_G_h, "(EPI_ISL_[0-9]+)$")[,1] ), 
               str_match( raw.id_N_h, "^[A-Z0-9]+")[,1] )
ac.raw_n <- c( gsub( "EPI_ISL_", "EPI", str_match( raw.id_G_n, "(EPI_ISL_[0-9]+)$")[,1] ),
               str_match( raw.id_N_n, "^[A-Z0-9]+")[,1] )

# match to id with tidy time format

raw_7326_ha <- fastaEx("data/out/pH5_7326.fasta")$id
raw_7138_na <- fastaEx("data/out/pNA_7138.fasta")$id

temp = c( "data/out//pH5_7326.fasta", "data/out//pNA_7138.fasta" )
for(i in 1:2)
{
  dir.i <- temp[i]
  id.i  <- fastaEx( dir.i )$id
  ac.i  <- str_match( id.i, "^[A-Z90-9]+" )[,1]
  yr.i  <- str_match( id.i, "_([0-9.]+)$" )[,2]
  
  if( grepl( "pH5", temp[i] ) )
  {
    m <- t.raw_h[ match( ac.i, ac.raw_h) ]
    if( TRUE %in% is.na(m) ){stop()}
    
  }else{
    
    m <- t.raw_n[ match( ac.i, ac.raw_n) ]
    if( TRUE %in% is.na(m) ){stop()}
  }
  
  nonmon   <- which( endsWith( m, "-99-99" ) )
  partialT <- rep( 0, length(m))
  partialT[nonmon] = 1
  
  df     <- data.frame( id = id.i, ac = ac.i, yr = yr.i, yr0 = m, partialT )
  write.table( df, 
               file  = paste0("./others/out/", i), 
               quote = FALSE, 
               row.names = FALSE, col.names = TRUE)
  
}  
