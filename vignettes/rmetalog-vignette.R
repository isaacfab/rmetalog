## ------------------------------------------------------------------------
library(rmetalog)

data("fishSize")
summary(fishSize)

## ----eval=FALSE----------------------------------------------------------
#  metalog()

## ------------------------------------------------------------------------
my_metalog <- metalog(fishSize$FishSize,
                      term_limit = 9,
                      bounds=0,
                      boundedness = 'sl',
                      step_len = .01)

## ------------------------------------------------------------------------
summary(my_metalog)

## ----message=FALSE, warning=FALSE, fig.width=6---------------------------
plot(my_metalog)

## ----message=FALSE, warning=FALSE, fig.width=6---------------------------
s <- rmetalog(my_metalog, n = 1000, term = 9)
hist(s)

## ----message=FALSE, warning=FALSE----------------------------------------
qmetalog(my_metalog, y = c(0.25, 0.5, 0.75), term = 9)

## ----message=FALSE, warning=FALSE----------------------------------------
pmetalog(my_metalog, q = c(3,10,25), term = 9)

## ----message=FALSE, warning=FALSE----------------------------------------
dmetalog(my_metalog, q = c(3,10,25), term = 9)

