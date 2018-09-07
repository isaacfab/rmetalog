## ----eval=FALSE----------------------------------------------------------
#  library(devtools)
#  install_github('isaacfab/RMetalog')
#  

## ------------------------------------------------------------------------
library(RMetalog)
data("fishSize")
summary(fishSize)

## ----eval=FALSE----------------------------------------------------------
#  r_metalog()

## ------------------------------------------------------------------------
my_metalog <- r_metalog(fishSize$FishSize,
                       term_limit = 9,
                       bounds=0,
                       boundedness = 'sl',
                       step_len = .01)

## ------------------------------------------------------------------------
r_metalog_summary(my_metalog)

## ----message=FALSE, warning=FALSE, fig.width=6---------------------------
r_metalog_plot(my_metalog)

## ----message=FALSE, warning=FALSE, fig.width=6---------------------------
s<-r_metalog_sample(my_metalog,n=1000,term=9)
hist(s)

## ----message=FALSE, warning=FALSE----------------------------------------
r_metalog_quantile(my_metalog,y=c(0.25,0.5,0.75),term=9)

