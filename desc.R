# In this script we extract descriptive statistics of the sample used for analysis of neuropsychological outcomes of cognitive
# rehabilitation including exporting tables for reporting in a final manuscript

# list packages to be used
pkgs <- c( "rstudioapi", "tidyverse", "dplyr", "openxlsx" )

# load or install each of the packages as needed
for (i in pkgs) {
  if (i %in% rownames (installed.packages()) == F) install.packages(i) # install if it ain't installed yet
  if ( i %in% names (sessionInfo() $otherPkgs) == F ) library(i, character.only = T ) # load if it ain't loaded yet
}

# set working directory (works in RStudio only)
setwd( dirname(getSourceEditorContext()$path) )

# create folders for figures, tables and session info
# prints TRUE and creates the folder if it was not present, prints NULL if the folder was already present.
sapply( c("figs", "tabs", "sess", "_nogithub", "_nogithub/data"), function(i) if( !dir.exists(i) ) dir.create(i) )

# read the data
d <- read.csv( "_nogithub/data/df.csv", sep = "," ) # the data set
rows <- c( "Age", "Years_of_education", read.csv( "_nogithub/data/liks.csv", sep = "," )$out ) # variables to describe


# ---- descriptive statistics mining ----

# prepare a data frame for descriptive stats
desc <- lapply( data.frame( n = "n", mean = "mean", sd = "sd", median = "median", min = "min", max = "max" ),
                # for each stat set-up a data.frame
                function(i)
                  matrix(
                    # set-up an empty matrix for each measure/stat combination
                    data = NA, nrow = length(rows), ncol = 6,
                    # give names to rows and columns
                    dimnames = list(
                      x = rows,
                      y = expand.grid( unique(d$assessment), unique(d$GROUP_EXP_CON) ) %>%
                        mutate( V3 = paste( Var2, Var1, sep = "_" ) ) %>%
                        select(V3) %>% unlist()
                    )
                  ) %>% as.data.frame() # re-format to a data.frame for easier manupulation down the line
)

# fill-in frequencies of the continuous variables
for ( i in names(desc$n) ) {
  
  grp <- strsplit(i,"_")[[1]] # prepare a grouping variable based on the column name
  
  # count non-zero observations for each variable/group pair (represented as n_females/n_males)
  desc$n[ rows , i ] <- paste0(
    apply( d[ with( d, GROUP_EXP_CON == grp[1] & assessment == grp[2] & Gender == "žena" ), rows ], 2, function(x) sum( !is.na(x) ) ), "/",
    apply( d[ with( d, GROUP_EXP_CON == grp[1] & assessment == grp[2] & Gender == "muž" ), rows ], 2, function(x) sum( !is.na(x) ) )
  )
  
  # calculate summary statistics
  # variables stay the same as for the counts above
  for ( j in names(desc)[-1] ) desc[[j]][ rows, i ] <- apply( d[ with( d, GROUP_EXP_CON == grp[1] & assessment == grp[2] ), rows ], 2, j, na.rm = T )
  
}

# save the descriptives as .xlsx
write.xlsx( desc, "tabs/desc.xlsx", rowNames = T )

# prepare a descriptive part of finalised tab for the text
t <- sapply( colnames(desc$n),
             # loop through all group/occasion combinations
             function(i)
               # print results as M (SD)
               paste0( sprintf( "%.2f", round( desc$mean[,i], 2 ) ), " (", sprintf( "%.2f", round( desc$sd[,i], 2 ) ), ")" )
             ) %>%
  # add variable names
  `rownames<-`( rows ) %>% as.data.frame() %>% rownames_to_column( "outcome" )

# save as .csv
write.table( t, "tabs/textab.csv", sep = ",", row.names = F, quote = F )


# ---- session info ----

# write the sessionInfo() into a .txt file
capture.output( sessionInfo(), file = "sess/desc.txt" )

