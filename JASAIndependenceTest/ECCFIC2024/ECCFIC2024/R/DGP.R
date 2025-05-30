#' @title Data generate process
#'
#' @param example example 1-10
#' @param model choose model, for example 1 and 7
#' @param n sample size
#' @param p dimension of \eqn{X}
#' @param q dimension of \eqn{Y}
#' @param delta parameter \eqn{\delta} in example 2
#' @param slice.num number of slices in example 4
#'
#' @return a list contain \eqn{X, Y}, and the discrete \eqn{Y}. 
#' @export
#'
#' @examples example = 1
#' model = 1
#' n = 25
#' p = 5
#' q = 1
#' dat = DGP(exmaple = example, model = model, n = 25, p = 5, q = 1)
DGP <- function(example, model = NULL, n = NULL, p = NULL, q = NULL, delta = NULL, slice.num = NULL){
  ## Use function switch() to choose example
  switch(
    example,
    {
      ## Example 1, Model (a)-(d)
      switch(
        model,
        {
          x = matrix(rnorm(n * p), n, p)
          y = matrix(rnorm(n * q), n, q)
        },
        {
          x = matrix(rt(n * p, df = 1), n, p)
          y = matrix(rt(n * q, df = 1), n, q)
        },
        {
          x = matrix(rt(n * p, df = 2), n, p)
          y = matrix(rt(n * q, df = 2), n, q)
        },
        {
          x = matrix(rt(n * p, df = 3), n, p)
          y = matrix(rt(n * q, df = 3), n, q)
        }
      )
      ## Quantile slicing method for making y discrete
      y.d <- cut(
        y, breaks = quantile(y, probs = seq(0, 1, by = 1/5)),
        labels = F, include.lowest = TRUE
      )
      return(list(x = as.matrix(x), y = as.matrix(y), y.d = as.matrix(y.d)))
    },
    ## Example 2, Different group indicator
    {
      x = replicate(p, c(
        rt(90, df = 4), rt(15, df = 4, ncp = delta),
        rt(15, df = 4, ncp = -delta)
      ))
      switch(
        model,
        {
          y.d = rep(1:4, each = 30)
        },
        {
          y.d = rep(c(1, 8, 0.2, 2.5), each = 30)
        }
      )
      return(list(x = as.matrix(x), y.d = as.matrix(y.d)))
    },
    ## Example 3, Unfinished
    {
      print('Unfinished')
    },
    ## Example 4, Aircraft data set in [sm] package
    {
      air.dat = sm::aircraft
      x = log(air.dat$Speed)
      y = log(air.dat$Span)
      ## Quantile slicing method for making y discrete
      y.d = cut(
        y, breaks = unique(quantile(y, probs = seq(0, 1, by = 1/slice.num))),
        labels = F, include.lowest = TRUE
      )
      return(list(x = as.matrix(x), y.d = as.matrix(y.d)))
    },
    ## Example 5, univariate continuous
    {
      x = rnorm(n)
      y = dnorm(x)
      slice.num = n / 5 * (n <= 20) + 5 * (n > 20)
      y.d = cut(
        y, breaks = quantile(y, probs = seq(0, 1, by = 1/slice.num)),
        labels = F, include.lowest = TRUE
      )
      return(list(x = as.matrix(x), y = as.matrix(y), y.d = as.matrix(y.d)))
    },
    ## Example 6, multivariate continuous
    {
      x = matrix(rnorm(n*p), n, p)
      y = x*matrix(rnorm(n*p), n, p)
      return(list(x = as.matrix(x), y = as.matrix(y)))
    },
    ## Example 7, regression setups A(1)-B(3)
    {
      beta = c(rep(1, 5), rep(0, 5))
      switch(
        model,
        ## A(1)
        {
          eps = rnorm(n)
          x = matrix(rnorm(n*10), n, 10)
          y = (x%*%beta)^2 + eps
          ## Quantile slicing method for making y discrete
          y.d <- cut(
            y, breaks = quantile(y, probs = seq(0, 1, by = 1/5)),
            labels = F, include.lowest = TRUE
          )
        },
        ## A(2)
        {
          eps = rnorm(n)
          x = cbind(
            rnorm(n, -8, 2), rf(n, 4, 10), rchisq(n, 5),
            rt(n, 15), rt(n, 3), matrix(rnorm(n*5), n, 5)
          )
          y = (x%*%beta)^2 + eps
          ## Quantile slicing method for making y discrete
          y.d <- cut(
            y, breaks = quantile(y, probs = seq(0, 1, by = 1/5)),
            labels = F, include.lowest = TRUE
          )
        },
        ## A(3)
        {
          eps = rnorm(n)
          x = cbind(
            matrix(rpois(n*5, 1), n, 5), matrix(rnorm(n*5), n, 5)
          )
          y = (x%*%beta)^2 + eps
          ## Quantile slicing method for making y discrete
          y.d <- cut(
            y, breaks = quantile(y, probs = seq(0, 1, by = 1/5)),
            labels = F, include.lowest = TRUE
          )
        },
        ## B(1)
        {
          eps = rnorm(n)
          x = matrix(rnorm(n*10), n, 10)
          y = 0.2*(x%*%beta)^2*eps
          ## Quantile slicing method for making y discrete
          y.d <- cut(
            y, breaks = quantile(y, probs = seq(0, 1, by = 1/5)),
            labels = F, include.lowest = TRUE
          )
        },
        ## B(2)
        {
          eps = rnorm(n)
          x = cbind(
            rnorm(n, -8, 2), rt(n, 5), rgamma(n, 9, 0.5),
            rf(n, 5, 12), rchisq(n, 5), matrix(rnorm(n*5), n, 5)
          )
          y = 0.2*(x%*%beta)^2*eps
          ## Quantile slicing method for making y discrete
          y.d <- cut(
            y, breaks = quantile(y, probs = seq(0, 1, by = 1/5)),
            labels = F, include.lowest = TRUE
          )
        },
        ## B(3)
        {
          eps = rnorm(n)
          x = cbind(
            matrix(rpois(n*5, 1), n, 5), matrix(rnorm(n*5), n, 5)
          )
          y = 0.2*(x%*%beta)^2*eps
          ## Quantile slicing method for making y discrete
          y.d <- cut(
            y, breaks = quantile(y, probs = seq(0, 1, by = 1/5)),
            labels = F, include.lowest = TRUE
          )
        }
      )
      return(list(x = as.matrix(x), y = as.matrix(y), y.d = as.matrix(y.d)))
    },
    ## Example 8
    {
      print('Unfinished')
    },
    ## Example 9
    {
      print('Unfinished')
    },
    ## Example 10
    {
      print('Unfinished')
    }
  )

}



