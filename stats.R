# In this script we fit series of models and summarise results. The modelling strategy is that of fitting GLMMs with neuropsychology
# outcome being predicted by group (experimental vs control), occassion (pre-test, post-test, follow-up) and their interaction

# list packages to be used
pkgs <- c( "rstudioapi", "tidyverse", "dplyr", "lme4", "lmerTest", "emmeans", "openxlsx" ) 

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

# read the data
d <- read.csv( "_nogithub/data/df.csv", sep = "," ) # the data set
l <- read.csv( "_nogithub/data/liks.csv", sep = "," ) # likelihoods to use
t <- read.csv( "tabs/textab.csv", sep = "," ) # descriptive part of the primary table


# ---- 2do ----

# exclude MDRS (subscales) from the analyses
# prepare big descriptive table in format: each cell = "mean (sd)", columns = "EXP1,EXP2,EXP3,CON1,CON2,CON3|stat.results (only p < .05)"
# add full-blown stat table for appendix



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
                rename_with( ~ "p value", contains("Pr(>|") )
              
              )

# finish the table
t2 <- do.call( rbind.data.frame, t2 ) %>%
  rownames_to_column( var = "Outcome" ) %>%
  mutate( Outcome = substr( Outcome, 1, nchar(Outcome)-2 ) ) %>%
  mutate( Transformation = sapply( Outcome, function(i) l[ l$out == i, "trans" ] ), .after = Outcome ) %>%
  mutate( Likelihood = sapply( Outcome, function(i) l[ l$out == i, "lik" ] ), .after = Outcome )


# ---- calculate contrasts ----

# extract 1) baseline between-group comparisons, 2) group/occasion interactions and 3) occasion differences per group simple effect
t3 <- lapply( setNames( c("b", "s", "i" ), c("baseline","simple","interaction") ),
              # loop through all contrasts
              function(i)
                lapply( setNames(l$out,l$out),
                        # loop through all outcomes
                        function(j) {
                          # fill-in baseline differences
                          if ( i == "b" ) return( pairs( emmeans( m[[j]], ~ GROUP_EXP_CON | assessment ), simple = "GROUP_EXP_CON" ) %>% as.data.frame() %>% `colnames<-` ( sub( "t.ratio|z.ratio", "test.stat", colnames(.) ) ) %>% filter( assessment == "1baseline" ) )
                          if ( i == "s" ) return( pairs( emmeans( m[[j]], ~ GROUP_EXP_CON | assessment ), simple = "assessment" ) %>% as.data.frame() %>% `colnames<-` ( sub( "t.ratio|z.ratio", "test.stat", colnames(.) ) ) )
                          if ( i == "i" ) return( contrast( emmeans( m[[j]], ~ GROUP_EXP_CON * assessment ), interaction = c("consec","revpairwise") ) %>% as.data.frame() %>% `colnames<-` ( sub( "t.ratio|z.ratio", "test.stat", colnames(.) ) ) )
                        }
                      ) %>%
                # pool all contrast per outcome
                do.call( rbind.data.frame, . ) %>% rownames_to_column( "outcome" ) %>%
                mutate( outcome = sub( "\\..*", "", outcome ), .before = 1 ) %>%
                mutate(`sig. (PCER = 5%)` =  ifelse( p.value > .05, NA, "*" ) )
              )

# extract Benjamini-Hochberg corrected threshold for a 5% FDR in interaction contrasts
bh_thres <- data.frame( p = sort( t3$interaction[ , "p.value"] ), # order the p-values from lowest to largest
                        thres = .05 * (1:nrow(t3$interaction) ) / nrow(t3$interaction) # prepare BH thresholds for each p-value
                        ) %>%
  # flag BH-significant p-values and extract the largest threshold (https://doi.org/10.1111/j.2517-6161.1995.tb02031.x)
  mutate( sig = ifelse( p <= thres, T, F ) ) %>% filter( sig == T ) %>% select(thres) %>% max()

# flag results that are statistically significant on 5% FDR
t3$interaction <- t3$interaction %>% mutate(`sig. (FDR = 5%)` = ifelse( p.value > bh_thres, NA, "*" ) )


# ---- add significant contrasts to the main table ----

# extract names of outcomes that had a p-value < .05 in at least one parameter
sigs <- with( t3, outcome[ `sig. (PCER = 5%)` == ":-)" ] ) %>% na.omit() %>% c() %>% unique()




# write the results into xlsx
write.xlsx( t2, "tabs/stats.xlsx", rowNames = T )


# ---- extract figures ----

# extract names of outcomes that had a p-value < .05 in at least one parameter
sigs <- with( t2, Outcome[ `sig. (PCER = 5%)` == ":-)" ] ) %>% na.omit() %>% c() %>% unique()

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
capture.output( sessionInfo(), file = "sess/pd_cogREH_stats_sessinfo.txt" )
