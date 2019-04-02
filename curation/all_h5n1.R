source("./function.R")

ls.gsgd        <- tagExtra( "./gsgd/out//pH5_7326_e0823.tre" )
ls.gsgd.wo2344 <- ls.gsgd$id[ which( ls.gsgd$tag == "00ff00" ) ]

leafEx( "./data/out/trim_pH5_7326.fasta", ls.gsgd.wo2344 )
rmDup( "curation/out/all_h5n1/trim_pH5_7326_4543.fasta", rmdup = TRUE)
rmdup_plus( "./curation/out/all_h5n1/trim_pH5_7326_4543_cr.fasta")

# removing confusing china isolates

reIntro.id  <- treeWalk( trefile = "./curation/out/all_h5n1/fasttree_pH5_7326_3_4543_cr_rmdP.tre", pattern = "00", showTree = TRUE)

rmDup( "./curation/out/all_h5n1/trim_pH5_7326_4543_cr_rmdP.fasta", 
       geo  = c( "China", "Hong_Kong"), 
       sero = "H5N1" )

leafEx( "./curation/out/all_h5n1/trim_pH5_7326_4543_cr_rmdP.fasta", 
        setdiff( fastaEx( "./curation/out/all_h5n1/pH5_7326_808.fasta" )$id, reIntro.id ) )


# remove outlier facilitated  with TempEst

outs <- c( "FJ602825_bar_headed_goose_Qinghai_10_2008_|China|_H5N1_2008.530",
           "FJ455820_environment_Qinghai_1_2008_|China|_H5N1_2008.445",
           "KX364460_swine_Shandong_SD1_2014_|China|_H5N1_2014.784",
           "FJ602820_bar_headed_goose_Qinghai_5_2008_|China|_H5N1_2008.456",
           "FJ390061_plateau_pika_Qinghai_04_2007_|China|_H5N1_2007.296",
           "KX364468_swine_Shandong_SD2_2014_|China|_H5N1_2014.838",
           "FJ602818_bar_headed_goose_Qinghai_3_2008_|China|_H5N1_2008.445" )

leafEx( "./curation/out/all_h5n1/pH5_7326_737.fasta", 
        setdiff( fastaEx( "./curation/out/all_h5n1/pH5_7326_737.fasta" )$id, outs ) ) 

cladeSampling( trefile = "./curation/out/all_h5n1/raxml_pH5_730e.tre", 
               fasfile = "./curation/out/all_h5n1/pH5_7326_730.fasta", 
               saveFasta = TRUE, suppList = FALSE, grid = 0.5)


# timeDice(fas.dir = "all_h5n1/pH5_730_cs.fasta", ecotable = FALSE, timetab.dir = "./time/time_ha_7326")
