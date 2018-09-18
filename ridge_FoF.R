ridgeFoF <- function (formula, data, subset, na.action, lambda = 0, model = FALSE, 
          x = FALSE, y = FALSE, contrasts = NULL,
          v_p = NULL, ...) 
{
  m <- match.call(expand.dots = FALSE)
  m$model <- m$x <- m$y <- m$contrasts <- m$... <- m$lambda <- m$v_p <- NULL
  m[[1L]] <- quote(stats::model.frame)
  m <- eval.parent(m)
  Terms <- attr(m, "terms")
  Y <- model.response(m)
  X <- model.matrix(Terms, m, contrasts)
  n <- nrow(X)
  p <- ncol(X)
  offset <- model.offset(m)
  if (!is.null(offset)) 
    Y <- Y - offset
  if (Inter <- attr(Terms, "intercept")) {
    Xm <- colMeans(X[, -Inter])
    Ym <- mean(Y)
    p <- p - 1
    X <- X[, -Inter] - rep(Xm, rep(n, p))
    Y <- Y - Ym
  }
  else{
    Xm <- colMeans(X)
    Ym <- mean(Y)
    X <- X - rep(Xm, rep(n, p))
    Y <- Y - Ym
  } 
  #Ym <- Xm <- NA
  Xscale <- drop(rep(1/(n), n) %*% X^2)^0.5 #n
  X <- X/rep(Xscale, rep(n, p))
  Xs <- svd(X)
  rhs <- t(Xs$u) %*% Y
  d <- Xs$d
  lscoef <- Xs$v %*% (rhs/d)
  lsfit <- X %*% lscoef
  resid <- Y - lsfit
  s2 <- sum(resid^2)/(n - p - Inter)
  HKB <- (p - 2) * s2/sum(lscoef^2)
  LW <- (p - 2) * s2 * n/sum(lsfit^2)
  if(is.null(v_p)){
    k <- length(lambda)
    dx <- length(d)
    div <- d^2 + rep(lambda, rep(dx, k))
    a <- drop(d * rhs)/div
    dim(a) <- c(dx, k)
    coef <- Xs$v %*% a
    dimnames(coef) <- list(names(Xscale), format(lambda))
    df <- colSums(matrix(d^2/div, dx))
  }else{
    v_lambda <- NULL
    for(i in 1:length(lambda)){
      v_lambda = append(v_lambda, rep(lambda[i], v_p[i]))
    }
    A = diag(d^2) + t(Xs$v)%*%diag(v_lambda)%*%(Xs$v)
    As = svd(A)
    a <- (As$v) %*% diag(1/(As$d)) %*% t(As$u) %*% diag(d) %*% rhs
    coef <- Xs$v %*% a
    dimnames(coef) <- list(names(Xscale))
    df <- sum(diag(diag(d^2) %*% (As$v) %*% diag(1/(As$d)) %*% t(As$u)))
  }
  GCV <- colSums((Y - X %*% coef)^2)/(n - df)^2
  res <- list(coef = drop(coef), scales = Xscale, Inter = Inter, 
              lambda = lambda, ym = Ym, xm = Xm, GCV = GCV, kHKB = HKB, 
              kLW = LW)
  class(res) <- "ridgelm"
  res
}