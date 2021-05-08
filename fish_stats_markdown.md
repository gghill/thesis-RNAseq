Statistical Analysis of Fish Condition
================
Griffin Hill
5/7/2021

#### Verification of sample inter-comparability for downstream RNA-seq comparisons

This notebook seeks to verify the validity of comparisons of *Boreogadus
saida* samples from fjord and offshore locations around the Greenland
Sea.

## Setup

A few libraries and other functions are critical to the analysis.

``` r
library(dplyr)
library(ggplot2)
library(cowplot)
library(stats)
library(scales)
# theme to add to ggplots to make text elements more legible
ax_theme = theme(
  axis.title.x = element_text(size = 16),
  axis.text.x = element_text(size = 14),
  axis.text.y = element_text(size = 14),
  axis.title.y = element_text(size = 16))

# function to clean up p values
fixp <- function(x, dig=3){
  x <- as.data.frame(x)
  
  if(substr(names(x)[ncol(x)],1,2) != "Pr")
    warning("The name of the last column didn't start with Pr. This may indicate that p-values weren't in the last row, and thus, that this function is inappropriate.")
  x[,ncol(x)] <- round(x[,ncol(x)], dig)
  for(i in 1:nrow(x)){
    if(x[i,ncol(x)] == 0)
      x[i,ncol(x)] <- paste0("< .", paste0(rep(0,dig-1), collapse=""), "1")
  }
  
  x
}

# function to calculate standard error
stderr <- function(x, na.rm=T) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}
```

## Import and mutation

The base dataframe includes fish identifiers, collection point, length
and weight, as well as a calculated Fulton K condition index.

``` r
all_fish = read.delim('./allfish_vis_ID.txt', header = T, dec = ',')
names(all_fish) = c("Fish_ID","GH_ID","location","length","weight","FultonK")
# more human legible names
all_fish[all_fish=="Bessel_Shelf"] <- "Bessel Offshore" 
all_fish[all_fish=="Shelf-N"] <- "North Shelf (Offshore)"
all_fish[all_fish=="Hoch"] <- "Hochsetter (Offshore)"
all_fish[all_fish=="Tyroler"] <- "Tyrolerfjord"
all_fish[all_fish=="Svalbard, Isfjord"] <- "Isfjord"
```

This dataframe can then be modified to introduce some additional useful
parameters for downstream analyses such as “type” which describes which
fjord-offshore pairing a sample belongs to and “enviro” which describes
simply whether its from a fjord or offshore environment (regardless of
specific
pairing).

``` r
# mutate to add a column of "type" describing which fjord-offshore pair each location belongs to
all_fish = mutate(all_fish,
                  type=ifelse(grepl("Bess",location,fixed=T),"Bess",
                              ifelse(grepl("Tyr",location,fixed=T)|grepl("Hoch",location,fixed=T), "Tyr", "Isfj"))
)
# column with habitat type
all_fish = mutate(all_fish,
                  enviro=ifelse(grepl("Fj",GH_ID,fixed=T),"Fjord","Offshore")
)
head(all_fish)
```

    ##    Fish_ID      GH_ID        location length weight FultonK type   enviro
    ## 1 Bsa17242 BesOff1701 Bessel Offshore    171  33.10    0.66 Bess Offshore
    ## 2 Bsa17244 BesOff1702 Bessel Offshore    166  26.50    0.58 Bess Offshore
    ## 3 Bsa17245 BesOff1704 Bessel Offshore    140  16.10    0.59 Bess Offshore
    ## 4 Bsa17246 BesOff1705 Bessel Offshore    156  21.30    0.56 Bess Offshore
    ## 5 Bsa17247 BesOff1706 Bessel Offshore    220  64.50    0.61 Bess Offshore
    ## 6 Bsa17248 BesOff1707 Bessel Offshore    204  49.75    0.59 Bess Offshore

## Basic visualization

Throughout this analysis, functions are used to iteratively create a
specific type of plot for all location pairings or habitat types.

``` r
histo_overview <- function(loc){
  ggplot(data=subset(all_fish, location == loc),aes(x = length, y = after_stat(count))) + 
    geom_histogram(bins = 10) +
    scale_y_continuous(breaks = c(1,5,10)) +
    coord_cartesian(ylim = c(0, 10)) +
    scale_x_continuous(breaks = c(120,160,200), limits = c(100,220)) +
    labs(y="count", x="length (mm)") +
    ggtitle(paste0(loc)) +
    theme(legend.position = 'none') +
    ax_theme
}
ordered_locs = c("Besselfjord", "Bessel Offshore", "Tyrolerfjord", "Hochsetter (Offshore)", "Isfjord", "North Shelf (Offshore)")
histo_list <- lapply(ordered_locs, histo_overview)
plot_grid(plotlist = histo_list,nrow = 3)
```

<img src="fish_stats_markdown_files/figure-gfm/length histos-1.png" style="display: block; margin: auto;" />

This approach combined with `lapply` and `plot_grid` allows for the
creation of ready made figures with grids of plots. A custom theme
profile and specific ordering of the locations provides figures with
legible labels and conveniently organizes the columns of the resulting
figure by location type (fjord on the left, offshore on the right).

This kind of count based visualization can also be achieved on a single
figure via grouping, but is more difficult to
decipher.

<img src="fish_stats_markdown_files/figure-gfm/unnamed-chunk-1-1.png" style="display: block; margin: auto;" />

## Length vs. weight relationship

Now, it makes sense to explore the relationship between length and
weight in the sampled fish a bit more. This is our first inkling of
actual condition of the organism, which of course we’d like to be fairly
uniform across individuals used in the study.

``` r
no_log_fish <- f + geom_point(aes(color = location))  +
  geom_smooth(method='lm',formula=y ~ x) +
  theme(legend.position = "bottom") +
  labs(x="length (mm)", y="weight (g)") +
  ggtitle("length weight relationship") +
  ax_theme
no_log_fish
```

<img src="fish_stats_markdown_files/figure-gfm/raw length weight-1.png" style="display: block; margin: auto;" />
It is apparent that this may not be the best fit as the regression is
linear and the data appears to be slightly curved. We can assess the fit
using the residuals with the goal being no pattern in the residuals
vs. fitted plot and similar variance over the whole range of the data.

``` r
par(mfrow=c(2,2))
plot(lm(all_fish$weight ~ all_fish$length, all_fish))
```

<img src="fish_stats_markdown_files/figure-gfm/no log resid-1.png" style="display: block; margin: auto;" />
Obviously our residuals have some shape and there may be a better way to
fit a model to this data. Based on the assumption of homogeneity of
variance over the entire sampled population, we can justify a log10
transformation of the
data.

``` r
log_fish <- ggplot(all_fish,aes(x=log10(length),y=log10(weight),location)) + geom_point(aes(color=location)) +
  geom_smooth(method = 'lm',formula=y~x) +
  theme(legend.position = 'bottom') +
  labs(x="log10(length (mm))", y="log10(weight (g))") +
  ggtitle("log10 transformed length weight relationship") +
  ax_theme
log_fish
```

<img src="fish_stats_markdown_files/figure-gfm/log transform-1.png" style="display: block; margin: auto;" />
This looks like a much better fit, but we can verify with the residual
plot again.

``` r
par(mfrow=c(2,2))
plot(lm(log10(all_fish$weight) ~ log10(all_fish$length),all_fish))
```

<img src="fish_stats_markdown_files/figure-gfm/log residuals-1.png" style="display: block; margin: auto;" />
Residuals look good as well, so now we can use this log transformation
going forward and break down our visualization by location.

``` r
log_lw_overview <- function(loc){
  ggplot(data=subset(all_fish, location == loc),aes(x = log10(length), y = log10(weight))) + 
    geom_point() +
    geom_smooth(method = 'lm',formula = y~x) +
    scale_y_continuous(limits = c(.5,2)) + # can comment out to scale or unscale
    scale_x_continuous(limits = c(2,2.4)) + # can comment out to scale or unscale
    labs(x="log10(length (mm))", y="log10(weight (g))") +
    ggtitle(paste0(loc)) +
    ax_theme
}
big_picture_list <- lapply(ordered_locs, log_lw_overview)
plot_grid(plotlist = big_picture_list,nrow = 3)
```

<img src="fish_stats_markdown_files/figure-gfm/all logged-1.png" style="display: block; margin: auto;" />

Scaling the axes in this plot allows for a better sense of the range of
values each location covers. We can see that different sites host a
different range of lengths, but all seem to fit our regression line
pretty well. Another way of representing this linear regression is
through its data
table.

``` r
out<-lapply(unique(all_fish$location) , function(z) { # log linear models
  data.frame(location=z, 
             coeff = coef(lm(log10(weight)~log10(length), data=all_fish[all_fish$location==z,]))[2],
             confint_low = confint.lm(lm(log10(weight)~log10(length), data=all_fish[all_fish$location==z,]))[2,1],
             confint_hi = confint.lm(lm(log10(weight)~log10(length), data=all_fish[all_fish$location==z,]))[2,2],
             n = count(subset(all_fish,all_fish$location==z)),
             row.names=NULL)
})
slopes = do.call(rbind ,out)
slopes
```

    ##                 location    coeff confint_low confint_hi  n
    ## 1        Bessel Offshore 2.986493    2.805167   3.167819 26
    ## 2            Besselfjord 2.855054    2.578276   3.131832 27
    ## 3           Tyrolerfjord 3.146002    2.941906   3.350098 30
    ## 4  Hochsetter (Offshore) 2.903567    2.743509   3.063624 28
    ## 5                Isfjord 2.994526    2.510008   3.479044 16
    ## 6 North Shelf (Offshore) 2.772301    2.625874   2.918728 31

Here we can see the coefficient (slope) for each regression line as well
as a confidence interval and n for each location. Beyond slope and range
of the regression line, Fulton’s K condition factor is a useful
indicator of similarity of condition across sites. Here we will compare
Fulton’s K (in which a higher value represents a healthier fish) and
slope (for which 3 represents the idealized cubic relationship between a
2 dimensional factor like length and a 3 dimensional one like weight).

``` r
# plot slope
slope_vis <- ggplot(slopes,aes(x=location,coeff,ymin=confint_low,ymax=confint_hi,color=location)) +
  geom_pointrange() + 
  geom_text(aes(label = paste0('n = ',n),y=confint_hi),position = position_nudge(y=.1)) +
  labs(y="coefficient") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y = element_text(size = 16),
        axis.text.y=element_text(size=14)) +
  theme(legend.position = 'none') +
  ggtitle('log-linear regression slope')
# plot Fulton's K
b = ggplot(all_fish,aes(location,FultonK,type))
all_K = b + geom_boxplot(aes(group=location, fill = location), outlier.color = 'red', outlier.shape = 1) +
  ggtitle("Fulton's K condition") + 
  theme(legend.position = 'none') + 
  theme(axis.title.x = element_blank()) +
  ax_theme +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))
# plot both on a shared axis
plot_grid(slope_vis,all_K,nrow=2,rel_heights = c(1,1.25))
```

<img src="fish_stats_markdown_files/figure-gfm/K and slope-1.png" style="display: block; margin: auto;" />

Here we can see visually that the confidence intervals appear to largely
overlap, especially those belonging to a fjord-offshore pairing, based
on both parameters. This suggests the locations are not significantly
different and indeed comparable.

Now that we’re thoroughly satisfied with that linear model output, we
can create a more detailed table that captures a complete description of
the model and its
outputs.

``` r
# list of function outputs that generates useful linear regression parameters for each location.
lm_out<-lapply((ordered_locs), function(z) { # log linear models
  data.frame(location=z, 
             n = count(subset(all_fish,all_fish$location==z)),
             log_intercept = summary(lm(log10(weight)~log10(length), data=all_fish[all_fish$location==z,]))['coefficients'][[1]][[2]],
             slope = coef(lm(log10(weight)~log10(length), data=all_fish[all_fish$location==z,]))[2],
             std_error = summary(lm(log10(weight)~log10(length), data=all_fish[all_fish$location==z,]))['coefficients'][[1]][[4]],
             r_squared = summary(lm(log10(weight)~log10(length), data=all_fish[all_fish$location==z,]))['adj.r.squared'],
             f_value = summary(lm(log10(weight)~log10(length), data=all_fish[all_fish$location==z,]))['fstatistic'][[1]][[1]],
             Pr_fvalue = anova(lm(log10(weight)~log10(length), data=all_fish[all_fish$location==z,]))[1,5],
             row.names=NULL)
}
)

lm_table = do.call(rbind,lm_out)

# use fixp() function to make p values human legible
lm_table = fixp(lm_table)
lm_table
```

    ##                 location  n log_intercept    slope  std_error adj.r.squared
    ## 1            Besselfjord 27      2.855054 2.855054 0.13438854     0.9454172
    ## 2        Bessel Offshore 26      2.986493 2.986493 0.08785616     0.9788050
    ## 3           Tyrolerfjord 30      3.146002 3.146002 0.09963650     0.9717064
    ## 4  Hochsetter (Offshore) 28      2.903567 2.903567 0.07786668     0.9809385
    ## 5                Isfjord 16      2.994526 2.994526 0.22427556     0.9268072
    ## 6 North Shelf (Offshore) 31      2.772301 2.772301 0.07159432     0.9803719
    ##     f_value Pr_fvalue
    ## 1  451.3403    < .001
    ## 2 1155.5227    < .001
    ## 3  996.9678    < .001
    ## 4 1390.4655    < .001
    ## 5  178.2757    < .001
    ## 6 1499.4211    < .001

This table can be easily exported as a .csv for use in a written work
using a command such as `write.table()`

However, we’ve so far only looked at one approximation of the
individual’s condition at a time (slope or K), what happens if these
are combined?

``` r
# remove one row with NA values
all_fish = all_fish[!all_fish$GH_ID=="IsFj1707",]
# define mean K values
mean_K = aggregate(all_fish[,6], list(all_fish$location), mean)
# calculate standard errors on these mean K values
se_K = aggregate(all_fish[,6], list(all_fish$location), stderr)
names(se_K) = c("location", "se_K")
mean_K = cbind(mean_K,se_K$se_K)
names(mean_K) = c("location","mean_K","se_K")
row.names(mean_K) = mean_K$location
# create tables in matching order of mean and standard error
row.names(slopes) = slopes$location
order = row.names(slopes)
mean_K = mean_K %>%
  slice(match(order, location))

# output is a table of means and confidence intervals/standard errors in 2 dimensions
plot_slope_K = cbind(slopes,mean_K$mean_K,mean_K$se_K)
head(plot_slope_K,1)
```

    ##                        location    coeff confint_low confint_hi  n
    ## Bessel Offshore Bessel Offshore 2.986493    2.805167   3.167819 26
    ##                 mean_K$mean_K mean_K$se_K
    ## Bessel Offshore     0.6080769   0.0070606

``` r
names(plot_slope_K) = c("location", "slope", "slope_low", "slope_hi", "n", "mean_K", "se_K")

# assign names and plot
slope_K <- ggplot(plot_slope_K,aes(x=slope,y=mean_K,xmin=slope_low,xmax=slope_hi,ymin=mean_K-se_K,ymax=mean_K+se_K,color=location)) +
  geom_pointrange() +
  geom_errorbarh() +
  ggtitle(paste0("Comparison of modeled slope and average measured","\n","Fulton's K value by location"))

slope_K
```

<img src="fish_stats_markdown_files/figure-gfm/slope and K-1.png" style="display: block; margin: auto;" />
From this we can see that there appears to be a lot more variability in
terms of mean Fulton’s K, even if we rescaled the y axis, some of the
confidence intervals wouldn’t overlap. There’s a more quantitative
approach we can take than visual comparison at this stage.

## Analysis of Variance (ANOVA) of Fulton’s K condition factor

First, we need to fit an analysis of variance model to the linear models
for each location. Fortunately, R does this on its own using the `aov()`
function. ANOVA assumptions: *Independence of the observations. *Each
subject should belong to only one group. *No significant outliers in any
cell of the design. *Normality. the data for each design cell should be
approximately normally distributed. \*Homogeneity of variances.

Based on the trawl based sampling regime, we can assume the observations
are independent and that each belongs to only one group. While some
outliers have emerged for individual sites, they fall within the normal
range of the whole population across sites, so it is reasonable to
assume that they do not significantly influence the results. Using our
initial histograms we can assume that the collection at each location is
roughly normal. The residuals from the log transformed data allow us to
assume homogeneity of variance.

``` r
aov_K_loc = aov(FultonK ~ location, data = all_fish)
summary(aov_K_loc)
```

    ##              Df Sum Sq Mean Sq F value   Pr(>F)    
    ## location      5 0.3410 0.06821    18.1 8.03e-14 ***
    ## Residuals   140 0.5277 0.00377                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 11 observations deleted due to missingness

With that highly significant p-value, we can conclude that some of the
Fulton’s K condition factors by location are different, but which ones?
For that, we need a test of significant differences. Since we are doing
pairwise comparisons, Tukey’s test for Honest Significant Differences
(HSD) is a good fit.

``` r
tuk_K_loc = TukeyHSD(aov_K_loc)
head(tuk_K_loc)
```

    ## $location
    ##                                                     diff          lwr
    ## Besselfjord-Bessel Offshore                   0.04414530 -0.004599095
    ## Hochsetter (Offshore)-Bessel Offshore         0.08712308  0.037431394
    ## Isfjord-Bessel Offshore                       0.09725641  0.039736934
    ## North Shelf (Offshore)-Bessel Offshore        0.01007123 -0.038673169
    ## Tyrolerfjord-Bessel Offshore                  0.13192308  0.082720979
    ## Hochsetter (Offshore)-Besselfjord             0.04297778 -0.006260752
    ## Isfjord-Besselfjord                           0.05311111 -0.004017337
    ## North Shelf (Offshore)-Besselfjord           -0.03407407 -0.082356425
    ## Tyrolerfjord-Besselfjord                      0.08777778  0.039033384
    ## Isfjord-Hochsetter (Offshore)                 0.01013333 -0.047805488
    ## North Shelf (Offshore)-Hochsetter (Offshore) -0.07705185 -0.126290382
    ## Tyrolerfjord-Hochsetter (Offshore)            0.04480000 -0.004891683
    ## North Shelf (Offshore)-Isfjord               -0.08718519 -0.144313634
    ## Tyrolerfjord-Isfjord                          0.03466667 -0.022852809
    ## Tyrolerfjord-North Shelf (Offshore)           0.12185185  0.073107458
    ##                                                      upr        p adj
    ## Besselfjord-Bessel Offshore                   0.09288969 9.995657e-02
    ## Hochsetter (Offshore)-Bessel Offshore         0.13681476 1.845678e-05
    ## Isfjord-Bessel Offshore                       0.15477589 4.047485e-05
    ## North Shelf (Offshore)-Bessel Offshore        0.05881562 9.910896e-01
    ## Tyrolerfjord-Bessel Offshore                  0.18112517 2.583500e-11
    ## Hochsetter (Offshore)-Besselfjord             0.09221631 1.247237e-01
    ## Isfjord-Besselfjord                           0.11023956 8.444251e-02
    ## North Shelf (Offshore)-Besselfjord            0.01420828 3.257007e-01
    ## Tyrolerfjord-Besselfjord                      0.13652217 1.002598e-05
    ## Isfjord-Hochsetter (Offshore)                 0.06807215 9.959013e-01
    ## North Shelf (Offshore)-Hochsetter (Offshore) -0.02781332 1.856292e-04
    ## Tyrolerfjord-Hochsetter (Offshore)            0.09449168 1.028115e-01
    ## North Shelf (Offshore)-Isfjord               -0.03005674 2.913951e-04
    ## Tyrolerfjord-Isfjord                          0.09218614 5.069206e-01
    ## Tyrolerfjord-North Shelf (Offshore)           0.17059625 4.451879e-10

``` r
gg_tuk = as.data.frame(tuk_K_loc$location)

# pick out pairings of interest to this study
interest = c("Besselfjord-Bessel Offshore", "Tyrolerfjord-Hochsetter (Offshore", "North Shelf (Offshore)-Isfjord")
gg_tuk = gg_tuk[interest,]
gg_tuk$pair = rownames(gg_tuk)


# plot pairwise TukeyHSD comparisons and color by significance level
ggplot(gg_tuk, aes(colour=cut(`p adj`, c(0, 0.01, 0.05, 1),
                           label=c("p<0.01","p<0.05","Non-Sig")))) +
  geom_hline(yintercept=0, lty="11", colour="grey30") +
  geom_errorbar(aes(pair, ymin=lwr, ymax=upr), width=0.2) +
  geom_point(aes(pair, diff)) +
  labs(colour="",x="Location pairing", y="Difference in modeled Fulton K value") +
  ggtitle("Fulton K Tukey HSD by location") +
  ax_theme +
  theme(legend.text  = element_text(size = 14)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
```

<img src="fish_stats_markdown_files/figure-gfm/Tukey HSD-1.png" style="display: block; margin: auto;" />
