# this script is used to summarize the sequence data used in this study

require( ggtree )


gisaid         = read.csv( "./data/in/G_2135_20171124.csv", stringsAsFactors = FALSE )
gisaid_isolate = gsub( "[A-Z]|_", "",  gisaid$Isolate_Id )
gisaid_HA      = str_extract( gisaid$HA.Segment_Id, "^[A-Z0-9]+" )

# Geo -------------------------------------------------------------------

tr1 <- read.nexus( "./beast/tree/geo_234_h5-anno.tre" )$tip.label
tr2 <- read.nexus( "./beast/tree/geo_232_h5-anno.tre" )$tip.label
tr3 <- read.nexus( "./beast/tree/geo_7_h5-anno.tre" )$tip.label
tr  <- c( tr1, tr2, tr3 )

tb1 <- read.table( "./beast/geo.234", header = TRUE )
tb2 <- read.table( "./beast/geo.232", header = TRUE )
tb3 <- read.table( "./beast/geo.7", header = TRUE )
tb  <- rbind( tb1, tb2, tb3 ) 

ac <- str_extract( tr, "^[A-Z0-9]+" )
ac[ grep( "EPI", ac ) ] = gisaid_HA[ match( gsub( "[A-Z]", "", grep( "EPI", ac, value = TRUE ) ), gisaid_isolate ) ]

name <- str_match( tr, "_([A-Za-z0-9_-]+)" )[,2]
name <- gsub( "_$", "", name )

state <- tb$states[ match(  tr, tb$id ) ]
state <- gsub( "^g", "", state )


df1 <- data.frame( accession_no  = ac,
                   sequence_name = name,
                   gene          = rep( "HA", length(ac) ),
                   sublineage    = c( rep( "2.3.4", length(tr1)), 
                                      rep( "2.3.2", length(tr2)), 
                                      rep( "7", length(tr3))),
                   state         = state,
                   analyses      = "Geo" )   


# All_H5N1 -------------------------------------------------------------------

tr <- read.nexus( "./beast/tree/all_h5n1-anno.tre" )$tip.label
ac <- str_extract( tr, "^[A-Z0-9]+" )
ac[ grep( "EPI", ac ) ] = gisaid_HA[ match( gsub( "[A-Z]", "", grep( "EPI", ac, value = TRUE ) ), gisaid_isolate ) ]

name <- str_match( tr, "_([A-Za-z0-9_-]+)" )[,2]
name <- gsub( "_$", "", name )


df2 <- data.frame( accession_no  = ac,
                   sequence_name = name,
                   gene          = rep( "HA", length(ac) ),
                   sublineage    = ".",
                   state         = ".",
                   analyses      = "All_H5N1" )  


# Sub_H5N1 -------------------------------------------------------------------

tr11 = read.tree( "./beast/dynamics/nwk/grid_232_h5.nwk" )$tip.label
tr12 = read.tree( "./beast/dynamics/nwk/grid_232_n1.nwk" )$tip.label
tr21 = read.tree( "./beast/dynamics/nwk/grid_234_h5.nwk" )$tip.label
tr22 = read.tree( "./beast/dynamics/nwk/grid_234_n1.nwk" )$tip.label
tr31 = read.tree( "./beast/dynamics/nwk/grid_7_h5.nwk" )$tip.label
tr32 = read.tree( "./beast/dynamics/nwk/grid_7_n1.nwk" )$tip.label
tr = c( tr11, tr12, tr21, tr22, tr31, tr32 )

ac <- str_extract( tr, "^[A-Z0-9]+" )
ac[ grep( "EPI", ac ) ] = gisaid_HA[ match( gsub( "[A-Z]", "", grep( "EPI", ac, value = TRUE ) ), gisaid_isolate ) ]

name <- str_match( tr, "_([A-Za-z0-9_-]+)" )[,2]
name <- gsub( "_$", "", name )

df3 <- data.frame( accession_no  = ac,
                   sequence_name = name,
                   gene          = c( rep( "HA", length(tr11) ), rep( "NA", length(tr12) ), 
                                      rep( "HA", length(tr21) ), rep( "NA", length(tr22) ), 
                                      rep( "HA", length(tr31) ), rep( "NA", length(tr32) )),
                   sublineage    = c( rep( "2.3.2", length(tr11)+length(tr12) ), 
                                      rep( "2.3.4", length(tr21)+length(tr22) ), 
                                      rep( "7", length(tr31)+length(tr32) )),
                   state         = ".",
                   analyses      = "Sub_H5N1" )  

# Eco -------------------------------------------------------------------

tr1 = read.nexus( "./beast/tree/eco_232_h5-anno.tre" )$tip.label
tr2 = read.nexus( "./beast/tree/eco_234_h5-anno.tre" )$tip.label
tr3 = read.nexus( "./beast/tree/eco_7_h5-anno.tre" )$tip.label
tr  <- c( tr1, tr2, tr3 )

tb1 <- read.csv( "./state/out/eco_232.csv", header = TRUE )
tb2 <- read.csv( "./state/out/eco_234.csv", header = TRUE )
tb3 <- read.csv( "./state/out/eco_7.csv", header = TRUE )
tb  <- rbind( tb1, tb2, tb3 ) 


ac <- str_extract( tr, "^[A-Z0-9]+" )
ac[ grep( "EPI", ac ) ] = gisaid_HA[ match( gsub( "[A-Z]", "", grep( "EPI", ac, value = TRUE ) ), gisaid_isolate ) ]

name <- str_match( tr, "_([A-Za-z0-9_-]+)" )[,2]
name <- gsub( "_$", "", name )

state <- as.character( tb$states[ match(  tr, tb$name ) ] )
state <- ifelse( startsWith( state, "D"),  "Domestic", "Wild" )

df4 <- data.frame( accession_no  = ac,
                   sequence_name = name,
                   gene          = rep( "HA", length(ac) ),
                   sublineage    = c( rep( "2.3.4", length(tr1)), 
                                      rep( "2.3.2", length(tr2)), 
                                      rep( "7", length(tr3))),
                   state         = state,
                   analyses      = "Eco" )  



pooltb = rbind( df1, df2, df3, df4 )
write.csv( pooltb, "./others/ac_table.csv", row.names = FALSE  )
