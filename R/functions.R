##-----------------------------------------------------------------------------
#'
#' vadim
#' Stephane Le Vu
#' Started 11/01/2024
#'
##-----------------------------------------------------------------------------

##---- libs ----
{
  library(survival)
  library(rms)
  library(ggplot2)
}


##---- parms ----
parms <- list(ncase = 1000, kctl = 2,
              pn = 0.1, or = 2,
              mean_x_exp = 21,
              mean_x_unexp = 28)
outlb <-  c(0.01, .99)

##---- make_mcc ----
## simulate dataset with known risk factors.
# 1:k controls, binary exposure,
# continuous variable (delay) that depends on binary exposure (in cases)
# so that exposed cases tend to have shorter delays
make_mcc <- function(
    ncase = 10,
    kctl = 2, ## nb of ctl per case
    mean_x_exp = 21, ## mean delay for exposed case
    sd_x_exp = 0.3, ## sd of log
    sd_x_unexp = 0.3,
    mean_x_unexp = 28, ## mean delay for unexposed case and ctls
    pn = .1, ## baseline prob of binary exposure (among ctl)
    or = 2, ## OR of binary exposure ( P(exp|case) = or * pn )
    seed = 321
){
  ## generate vectors
  set.seed(seed)
  set <- y <-  b <- x <- numeric(0) ## initialize
  for (s in 1:ncase){ ## each stratum
    # case and ctl are exposed to b or not
    b_case <- rbinom(1, 1, or * pn)
    b_ctl <- rbinom(kctl, 1, pn)
    # continuous value for case and ctl
    x_case <- ifelse(b_case,
                     rlnorm(1, log(mean_x_exp), sd_x_exp ),
                     rlnorm(1, log(mean_x_unexp), sd_x_unexp ))
    x_ctl <- rlnorm(kctl, log(mean_x_unexp), sd_x_unexp )
    ## append values
    y <- c(y, 1, rep(0, kctl))
    set <- c(set, rep(s, 1 + kctl))
    b <- c(b, b_case, b_ctl)
    x <- c(x, round(x_case), round(x_ctl))
  }
  ## data frame
  df <- data.frame(set, y, b, x)
  df
}

##---- run sim ----
# df <- make_mcc(ncase = 1000, kctl = 5)
df <- do.call(make_mcc, parms)

##---- check ----
if (TRUE){
  head(df, 10)
  # raw OR
  m <- table(case = df$y, exposed = df$b)
  m <- matrix( as.numeric(m), ncol = ncol(m), dimnames = dimnames(m))
  or <- m[1,1] * m[2,2] / m[2,1] / m[1,2]
  ## plot
  df$cc <- ifelse(df$y == 1,'Case','Ctrl')
  df$expo <- ifelse(df$b == 1,'Exposed','Unexposed')
  p0 <- ggplot(data = df, aes(x, colour = expo, fill = expo)) +
    geom_density(alpha = .2) +
    facet_grid(vars(cc))
}

##---- MAR intervals ----
if (TRUE){
  fmar <- 0.2 # prop missing at random
  r <- sample(1:nrow(df), round(nrow(df) * fmar))
  df[r, 'x'] <- NA
}

##---- quantiles ----
{
  # discard interval outliers
  q <- quantile(df$x, probs = outlb, na.rm = TRUE)
  x <- with(df, ifelse(x > q[1] & x < q[2], x, NA))
  # categorize x
  nqt <- 5
  co <- quantile(x, probs = seq(0, 1, 1/nqt), na.rm = TRUE)
  df$xc <- cut(x, co, right = FALSE)
  # new expo variable
  df$expo <- factor(ifelse(df$b == 1,
                 paste0("E_", df$xc),
                 "UNEXP"),
                 levels = c("UNEXP", paste0("E_", levels(df$xc))) )
  ## Fit clogit on categories
  {
    m <- clogit(y ~ expo + strata(set), data = df)
    res <- as.data.frame(
      round( cbind(exp(coef(m)), exp(confint(m))), 2) )
    colnames(res) <- c("OR", "OR_lo", "OR_up")
  }

}

##---- spline ----
{
  # Fit conditional logistic regression with RCS and an interaction
  dd <- datadist(df)
  options(datadist='dd')
  sfit <- lrm(y ~ rcs(x, 4) * b, data = df, x = TRUE, y = TRUE)


  ## explicit prediction
  {
    bound <- round(quantile(df$x, probs = outlb, na.rm = TRUE)) ## rm 1% extreme
    xref <- parms$mean_x_unexp
    xs <- bound[1]:bound[2]
    a <- list(x = xs, b = 1)
    b <- list(x = xref, b = 0)
    con <- contrast(sfit, a, b)
    df1 <- data.frame(x = xs, lapply(con[c('Contrast','Lower','Upper')], exp), 1 )
    names(df1) <- c("x", "or", "or_lo", "or_up", "group")
  }


  ## plot
  p2 <- ggplot(df1, aes(x, or)) +
    geom_ribbon(aes(ymin = or_lo, ymax = or_up,
                    fill = "red"), alpha = 0.5)  +
    geom_hline(yintercept = 1, col = "grey") +
    geom_line(aes(colour = "red"), linewidth = 1) +
    # xlim( li[[1]] ) +
    xlab("Days since previous dose") +
    ylab("OR (95%CI)") +
    theme(legend.position = "none")
  p2
}
