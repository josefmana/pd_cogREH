# In this script we fit series of models and summarise results. The modelling strategy is that of fitting GLMMs with neuropsychology
# outcome being predicted by group (experimental vs control), occassion (pre-test, post-test, follow-up) and their interaction

# Run only after importing raw data via '00_dataprep.R' and calculating descriptives vi '01_desc.R'

# list packages to be used
pkgs <- c("rstudioapi","tidyverse","dplyr","purrr","lme4","lmerTest","emmeans","openxlsx","gt")

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

# set-up plotting theme
theme_set( theme_minimal( base_size = 18 ) )

# in-house function for printing rounded numbers
rprint <- function(x,dec) sprintf( paste0("%.",dec,"f"), round(x,dec) )

# read the data
d <- read.csv( "_nogithub/data/df.csv", sep = "," ) # the data set
l <- read.csv( "_nogithub/data/liks.csv", sep = "," ) # likelihoods to use
t <- read.csv( "tabs/textab.csv", sep = "," ) # descriptive part of the primary table


# ---- data pre-processing ----

# check that there's exactly one observation (or lack of thereof) in each patient/assessment pair
with( d, table( Code, assessment) ) %>% as.data.frame() %>% filter( Freq != 1 ) # print cells with counts different from one

# transform relevant variables:
# (i) log transform reaction times
# (ii) for binomial variables prepare a column with number of failures
# (iii) in-sample standardise Gaussian outcomes
for ( i in c( with( l, out[trans == "log"] ) %>% na.omit() ) ) d[ ,i] <- log( d[ ,i] )
for ( i in c( with( l, out[lik == "binomial"] ) %>% na.omit() ) ) d[ ,i] <- d[ ,i] / with( l, ceil[out==i] )
for ( i in c( with( l, out[lik == "gaussian"] %>% na.omit() ) ) ) d[ ,i] <- ( d[ ,i] - mean( d[ ,i], na.rm = T ) ) / sd( d[ ,i], na.rm = T )


# ---- models fitting ----

# prepare a set of formulas for fitting
f <- lapply( setNames( l$out, l$out ), function(i) paste0( i, " ~ 1 + assessment * GROUP_EXP_CON + (1 | Code)" ) %>% as.formula )

# let's fit it
m <- lapply( setNames( names(f), names(f) ),
             function(i) {
               if ( with( l, lik[ out == i ] ) == "gaussian" ) return( lmer( f[[i]], data = d, REML = F ) )
               if ( with( l, lik[ out == i ] ) == "poisson" ) return( glmer( f[[i]], data = d, family = poisson() ) )
               if ( with( l, lik[ out == i ] ) == "binomial" ) return( glmer( f[[i]], data = d, family = binomial(), weights = rep( with( l, ceil[out == i] ) , nrow(d) ) ) )
             }
          )


# --- model parameters ----

# extract parameter estimates from each model
t2 <- lapply( setNames( l$out, l$out ),
              
              # loop through outcomes to extract coefficient estimates and tidy up the outcome for binding the results
              function(i) summary(m[[i]])$coefficients %>%
                as.data.frame %>% # re-format so that the ensuing functions work
                rownames_to_column( var = "Predictor" ) %>%
                mutate( df = ifelse( with( l, lik[out == i] ) != "gaussian", NA, df ), .after = "Std. Error" ) %>% # add column for df for non-gaussian models
                rename_with( ~ "Test Stat", contains("value") ) %>% # unity column names for Z and t-tests
                rename_with( ~ "p.value", contains("Pr(>|") )
              
              ) %>%
  
  # finish the table
  do.call( rbind.data.frame, . ) %>%
  rownames_to_column( var = "Outcome" ) %>%
  mutate( Outcome = substr( Outcome, 1, nchar(Outcome)-2 ) ) %>%
  mutate( Transformation = sapply( Outcome, function(i) l[ l$out == i, "trans" ] ), .after = Outcome ) %>%
  mutate( Likelihood = sapply( Outcome, function(i) l[ l$out == i, "lik" ] ), .after = Outcome ) %>%
  mutate(`sig. (PCER = 5%)` =  ifelse( p.value > .05, NA, ":-)" ) )


# ---- calculate contrasts ----

# extract 1) simple between-group comparisons, 2) occasion per group simple comparisons and 3) group/occasion interactions
t3 <- lapply( setNames( c("g", "o", "i" ), c("simple_group","simple_occasion","interaction") ),
              # loop through all contrasts
              function(i)
                lapply( setNames(l$out,l$out),
                        # loop through all outcomes
                        function(j) {
                          # fill-in baseline differences
                          if ( i == "g" ) return( pairs( emmeans( m[[j]], ~ GROUP_EXP_CON | assessment ), simple = "GROUP_EXP_CON", reverse = T, type = "response" ) )
                          if ( i == "o" ) return( pairs( emmeans( m[[j]], ~ GROUP_EXP_CON | assessment ), simple = "assessment", reverse = T, type = "response" ) )
                          if ( i == "i" ) return( contrast( emmeans( m[[j]], ~ GROUP_EXP_CON * assessment ), interaction = c("consec","revpairwise"), type = "response" ) )
                        } 
                      )
              )

# loop through t3 and format all subtables
for ( i in names(t3) ) {
  
  # align column names withing each contrast for all outcomes
  for ( j in l$out ) t3[[i]][[j]] <- t3[[i]][[j]] %>%
      as.data.frame() %>%
      select( !contains("null") ) %>% # added for the ORs
      `colnames<-` ( c("contrast","condition","estimate","SE","df","test.stat","p.value") )
  
  # create a single table for each contrast
  t3[[i]] <- t3[[i]] %>%
    
    do.call( rbind.data.frame, . ) %>%
    rownames_to_column( "outcome" ) %>% # add outcome name
    mutate( outcome = sub( "\\..*", "", outcome ), .before = 1 ) %>% # format outcome name
    mutate( outcome = ifelse( outcome == "Total_", "Total_.AVLT.1_5", outcome ) ) %>% # return RAVLT total name back
    mutate(`sig. (PCER = 5%)` =  ifelse( p.value > .05, NA, ":-)" ) ) %>% # add per comparison error rate (PCER) 5% significance flags
    
    # add for text/appendix table information
    mutate(
      
      across( all_of( c("estimate","SE","df") ), ~ rprint(.x,2) ),
      test.stat = rprint(test.stat,3),
      p.value = ifelse( p.value < .001, " < .001", paste0( "= ", substr( rprint(p.value,3), 2, 5 ) ) ),
      
      contr = case_when(
        grepl( " - ", contrast ) ~ paste0( "b = ", estimate," (",SE,"), ","t(",df,") = ",test.stat,", p ",p.value ),
        grepl( " / ", contrast ) ~ paste0( "OR = ", estimate," (",SE,"), ","Z = ",test.stat,", p ",p.value )
      ) 
    )

}

# ---- keeping FDR correction out of the way as there was no significant results on 5% level (just report as unclear) ----

# extract Benjamini-Hochberg corrected threshold for a 5% FDR for interaction contrasts
#bh_thres <- data.frame( p = sort( t3$interaction[ , "p.value"] ), # order the p-values from lowest to largest
#                        thres = .05 * (1:nrow(t3$interaction) ) / nrow(t3$interaction) # prepare BH thresholds for each p-value
#                        ) %>%
#  # flag BH-significant p-values and extract the largest threshold (https://doi.org/10.1111/j.2517-6161.1995.tb02031.x)
#  mutate( sig = ifelse( p <= thres, T, F ) ) %>% filter( sig == T ) %>% select(thres) %>% max()
#
# flag results that are statistically significant on 5% FDR
#t3$interaction <- t3$interaction %>% mutate(`sig. (FDR = 5%)` = ifelse( p.value > bh_thres, "-", ":-)" ) )


# ---- appendix table with simple contrast statistics ----

# prepare a formatted table with all descriptive and inferential stats
t4 <- list(
  
  # stat table for simple occasion contrasts
  t3$simple_occasion %>%
    mutate(
      contrast =
        case_when(
          contrast %in% paste0("2retest ",c("-","/")," 1baseline") ~ "21",
          contrast %in% paste0("3followup ",c("-","/")," 2retest") ~ "31",
          contrast %in% paste0("3followup ",c("-","/")," 1baseline") ~ "32"
        )
    ) %>%
    select(outcome,condition,contrast,contr) %>%
    pivot_wider( values_from = contr, names_from = c(condition,contrast), names_prefix = "occasion_" ),
  
  # stat table for simple group contrasts
  t3$simple_group %>%
    select(outcome,condition,contr) %>%
    pivot_wider( values_from = contr, names_from = condition, names_prefix = "group_" )

) %>%
  
  reduce( full_join ) %>%
  left_join( select( t, 1,5:7, 2:4 ), ., by = "outcome" ) %>%
  
  # format it
  gt() %>%
  cols_align( columns = -1, align = "center" ) %>%
  
  tab_spanner( label = "EXP", columns = starts_with("EXP"), id = "EXPdesc", gather = F ) %>%
  tab_spanner( label = "CON", columns = starts_with("CON"), id = "CONdesc", gather = F ) %>%
  tab_spanner( label = "Descriptive statistics", columns = 2:7, gather = F ) %>%
  
  tab_spanner( label = "CON", columns = starts_with("occasion_CON"), id = "CONmod", gather = F ) %>%
  tab_spanner( label = "EXP", columns = starts_with("occasion_EXP"), id = "EXPmod", gather = F ) %>%
  tab_spanner( label = "Occasion simple contrasts", columns = starts_with("occasion") ) %>%
  
  tab_spanner( label = "Group simple contrasts", columns = starts_with("group") ) %>%
  tab_spanner( label = "Generalized Linear Mixed Models' results", columns = 8:16 ) %>%
  
  cols_label(
    ends_with("1baseline") ~ "Baseline",
    ends_with("2retest") ~ "Re-test",
    ends_with("3followup") ~ "Follow-up",
    ends_with("21") ~ "Re-test-minus-Baseline",
    ends_with("31") ~ "Follow-up-minus-Baseline",
    ends_with("32") ~ "Follow-up-minus-Re-test",
    "outcome" ~ ""
  ) %>%
  
  tab_footnote(
    locations = cells_column_spanners("Descriptive statistics"),
    footnote = "Data are presented as mean ± standard deviation for continuous variables and as median (interquantile range) for count variables."
  ) %>%
  tab_footnote(
    locations = cells_column_spanners("Generalized Linear Mixed Models' results"),
    footnote = "Results indicate expected difference (b) or odds ratio (OR) with resulting test statistic (t-statistic with degrees of freedom in brackets or Z statistic) and p-values."
  ) %>%
  tab_footnote( locations = cells_column_spanners("Group simple contrasts"), footnote = "EXP-minus-CON." )

# save it
gtsave( t4, "tabs/bigtab.docx" )
  

# ---- add significant contrasts to the main table ----

# extract texts for all PCER significant contrasts
txt <- full_join( do.call( rbind.data.frame, t3[1:2] ) %>% add_column( type = "simple" ),
                  t3$interaction %>% add_column( type = "interaction" )
                  ) %>%
  
  # add reports of significant effects
  mutate(
    txt = case_when(
      type == "interaction" ~ case_when(
        `p.value` < .05 & estimate < 0 ~ paste0( "[ (", contrast, ") * (", condition, " ) ] < 0" ),
        `p.value` < .05 & estimate > 0 ~ paste0( "[ (", contrast, ") * (", condition, " ) ] > 0" )
      ),
      type == "simple" ~ case_when(
        `p.value` < .05 & estimate < 0 ~ paste0( sub( "-", "<", contrast ), " | ", condition ),
        `p.value` < .05 & estimate > 0 ~ paste0( sub( "-", ">", contrast ), " | ", condition )
      )
    )
  ) %>%
  
  # keep only the texts for further procesing
  select( outcome, txt ) %>%
  filter( complete.cases( txt ) )

# add likelihoods, transformations and significant contrasts summary
t <- t %>%
  # for each column, loop through all outcomes but age and education which were not compared
  mutate( likelihood = sapply( outcome, function(i) with( l, lik[out == i] ) ),
          transformation = sapply( outcome, function(i) with( l, trans[out == i] ) ),
          sig_results = sapply( outcome, function(i) with( txt, paste( txt[outcome == i], collapse = "\n" ) ) )
          )


# ---- export tables ----

# write the results into xlsx
write.xlsx( list( texttab = t,
                  parameters = t2,
                  interactions = t3$interaction,
                  simple_group = t3$simple_group,
                  simple_occasion = t3$simple_occasion
                  ),
            
            "tabs/stats.xlsx", rowNames = F )


# ---- extract figures ----

# extract names of outcomes that had a p-value < .05 in at least one parameter
sigs <- with( full_join( t3$interaction, do.call( rbind.data.frame, t3[1:2] ) ), outcome[ `sig. (PCER = 5%)` == "*" ] ) %>% na.omit() %>% c() %>% unique()

# prepare folders for emmeans
sapply( c("emmip", "comp" ), function(i) if( !dir.exists( paste0("figs/",i) ) ) dir.create( paste0("figs/",i) ) )

# save a figure per significant outcome
for ( i in sigs ) {
  
  # plot emmeans as a typical interaction plot
  emmip( m[[i]], GROUP_EXP_CON ~ assessment, CIs = T, type = "response" ) +
    labs( title = i ) +
    theme( legend.position = "bottom", plot.title =  element_text( hjust = .5 ) )
  ggsave( paste0( "figs/emmip/",i,".jpg"), dpi = 300, width = 13.1/1.5, height = 7.24 )
  
  # plot emmeans with per group comparisons via CIs adjusted by Tukey
  plot( emmeans( m[[i]], ~ assessment * GROUP_EXP_CON ), comparisons = T, by = "GROUP_EXP_CON" ) +
    coord_flip() +
    labs( title = i ) +
    theme( plot.title =  element_text( hjust = .5 ) )
  ggsave( paste0( "figs/comp/",i,".jpg"), dpi = 300, width = 13.1/1.5, height = 7.24 )
  
}


# ---- session info ----

# write the sessionInfo() into a .txt file
capture.output( sessionInfo(), file = "sess/stats.txt" )

