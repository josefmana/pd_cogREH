# In this script we prepare data set of neuropsychological outcomes after cognitive rehabilitation in a suitable format.

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
d <- list(
  bas = read.xlsx( "_nogithub/raw/PN celkem_exp_con_final.xlsx", sheet = "VSTUP exp + con" ),
  ret = read.xlsx( "_nogithub/raw/PN celkem_exp_con_final.xlsx", sheet = "Retest po 3 měsících" ),
  fol = read.xlsx( "_nogithub/raw/PN celkem_exp_con_final.xlsx", sheet = "Retest po 1 roce" )
)

# read a list of included patients
inc <- read.xlsx( "_nogithub/raw/PN celkem_exp_con_final.xlsx", sheet = "VSTUP_RETEST PO 3 MESICICH" )$A_Code


# ---- raw-to-data transformation ----

# do some housekeeping with column names
for ( i in names(d) ) {
  names(d[[i]]) <- names(d[[i]]) %>% sub("^[^_]*_", "", . ) # get rid of the prefixes
  d[[i]][ , "assessment"] <- case_when( i == "bas" ~ "1baseline", i == "ret" ~ "2retest", i == "fol" ~ "3followup" ) # add a variable denoting the assessment type
}

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

# rename DRS-2 scores such that they don't include slashes
names(d0)[ grepl("MDRS",names(d0)) ] <- ( names(d0)[ grepl("MDRS",names(d0)) ] %>% strsplit( "_(?!.*_)", perl = T ) %>% do.call( cbind, . ) )[1, ]

# keep only included patients
d1 <- d0[ d0$Code %in% inc , ]

# write as .csv 
write.table( d1, "_nogithub/data/df.csv", sep = ",", row.names = F, quote = F )


# ---- session info ----

# write the sessionInfo() into a .txt file
capture.output( sessionInfo(), file = "sess/dataprep.txt" )
