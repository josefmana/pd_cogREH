# In this script we first prepare tables for a suitable format then pre-porocess the data, fit series of models and summarise
# the results. The modelling strategy is that of fitting LMMs with increasingly more complex structure of predictors with
# subsequent model comparison for evalauting the importance of added predictors above and beyond predictors already included.

# list packages to be used
pkgs <- c( "rstudioapi", # setting working directory via RStudio API
           "tidyverse", "dplyr", # data wrangling
           "lme4", # model fitting
           "performance", # model quality checking
           "openxlsx" # reading and writing .xlsx
           )

# load or install each of the packages as needed
for (i in pkgs) {
  if (i %in% rownames (installed.packages()) == F) install.packages(i) # install if it ain't installed yet
  if ( i %in% names (sessionInfo() $otherPkgs) == F ) library(i, character.only = T ) # load if it ain't loaded yet
}

# set working directory (works in RStudio only)
setwd( dirname(getSourceEditorContext()$path) )

# create folders for figures, tables and session info
# prints TRUE and creates the folder if it was not present, prints NULL if the folder was already present.
sapply( c("figs", "tabs", "sess"), function(i) if( !dir.exists(i) ) dir.create(i) )


# ---- pre-processing  ----

# read the data
d <- list(
  bas = read.xlsx( "_no_github/PN celkem_exp_con_final.xlsx", sheet = "VSTUP exp + con" ),
  ret = read.xlsx( "_no_github/PN celkem_exp_con_final.xlsx", sheet = "Retest po 3 měsících" ),
  fol = read.xlsx( "_no_github/PN celkem_exp_con_final.xlsx", sheet = "Retest po 1 roce" )
)

# get rid of the prefixes and add a variable denoting the assessment type
for ( i in names(d) ) {
  names(d[[i]]) <- names(d[[i]]) %>% sub("^[^_]*_", "", . )
  d[[i]][ , "assessment"] <- case_when( i == "bas" ~ "baseline", i == "ret" ~ "retest", i == "fol" ~ "follow-up" )
}

# because baseline TMT-B scores were read as character rather than numeric, re-format them
d$bas$TMT_B_s <- as.numeric( d$bas$TMT_B_s )

# check whether column names coincide across data sets
for ( i in 1:3 ) print( isTRUE( all.equal( names( combn( names(d),2)[1,i] ), names( combn( names(d),2)[2,i] ) ) ) )
# it's A-OK

# collapse the sub data sets into a single long-format data set
d0 <- do.call( rbind.data.frame, d ) %>%
  arrange( by = Code ) %>% `rownames<-`( 1:nrow(.) ) %>% # rearrange and rename rows
  select( -c(2:4,6,9) ) %>% # keep only variables of interest
  mutate( GROUP_EXP_CON = case_when( GROUP_EXP_CON == 1 ~ "EXP", GROUP_EXP_CON == 2 ~ "CON" ) ) %>%
  # put a ceiling to reaction times scores
  mutate( TMT_A_s = ifelse( TMT_A_s > 360, 360, TMT_A_s ),
          TMT_B_s = ifelse( TMT_B_s > 600, 600, TMT_B_s ),
          Stroop_points = ifelse( Stroop_points > 60, 60, Stroop_points ),
          Stroop_neutral_words = ifelse( Stroop_neutral_words > 60, 60, Stroop_neutral_words ),
          Stroop_color_words = ifelse( Stroop_color_words > 240, 240, Stroop_color_words )
          )


# ---- data set description ----

# prepare a data frame for descriptive stats
desc <- lapply( data.frame( n = "n", mean = "mean", sd = "sd", median = "median", min = "min", max = "max" ),
                # for each stat set-up a data.frame
                function(i)
                  matrix(
                    # set-up an empty matrix for each measure/stat combination
                    data = NA, nrow = length( names(d0)[-c(1,2,ncol(d0))] ), ncol = 6,
                    # give names to rows and columns
                    dimnames = list(
                      x = names(d0)[-c(1,2,ncol(d0))],
                      y = expand.grid( unique(d0$assessment), unique(d0$GROUP_EXP_CON) ) %>%
                        mutate( V3 = paste( Var2, Var1, sep = "_" ) ) %>%
                        select(V3) %>% unlist()
                    )
                  ) %>% as.data.frame() # re-format to a data.frame for easier manupulation down the line
                )

# start by filling-in Gender data (fill-in frequencies only)
for ( i in names(desc$n) ) {
  
  grp <- strsplit(i,"_")[[1]] # prepare a grouping variable based on the column name
  desc$n[ "Gender", i ] <- with( d0, table( Gender[ GROUP_EXP_CON == grp[1] & assessment == grp[2] ] ) ) %>% paste(., collapse = "/" ) # fill-in gender distributions

} 

# next fill-in frequencies of the continuous variables
for ( i in names(desc$n) ) {
  
  grp <- strsplit(i,"_")[[1]] # prepare a grouping variable based on the column name
  vrs <- rownames(desc$n)[-1] # variables to be described
  
  # count non-zero observations for each variable/group pair
  desc$n[ vrs , i ] <- apply( d0[ with( d0, GROUP_EXP_CON == grp[1] & assessment == grp[2] ), vrs ], 2, function(x) sum( !is.na(x) ) )
  
  # calculate summary statistics
  # variables stay the same as for the counts above
  for ( j in names(desc)[-1] ) desc[[j]][ vrs, i ] <- apply( d0[ with( d0, GROUP_EXP_CON == grp[1] & assessment == grp[2] ), vrs ], 2, j, na.rm = T )
  
}

# save the descriptives as .xlsx
write.xlsx( desc, "tabs/pd_cogREH_desc.xlsx", rowNames = T )


# ---- session info ----

# write the sessionInfo() into a .txt file
capture.output( sessionInfo(), file = "sess/pd_cogREH_stats_sessinfo.txt" )
