#' Title
#'
#' @param Y outcome variable vector
#' @param A missingness indicator vector
#' @param piX covariate matrix or vector in the propensity score model
#' @param phiX covariate matrix or vector in the outcome model
#' @param initEta initial value of tilde_eta
#' @param ctrl a list of control parameters for the fitting of the outcome
#'
#' @return betaDR
#' @export
#'
#' @examples
DR <- function( Y, A, piX, phiX, initEta=NULL, ctrl=list(epsilon=1e-8,maxit=50) ) {
  if ( !all( A %in% c(0,1) ) ) stop( "íAí should be binary." )
  YisBin <- all( Y %in% c(0,1) )
  famY <- ifelse( YisBin, "binomial", "gaussian" )
  piMLE <- fitted( glm( A ~ piX, family=binomial ) )
  AdPi <- A/piMLE
  AdPi[ A==0 ] <- 0
  YAdPi <- Y * AdPi
  YAdPi[ A==0 ] <- 0
  phiX1 <- cbind( 1, phiX )
  if ( is.vector( phiX ) ) nEta <- 2
  else nEta <- dim(phiX)[[2]]+1
  if ( is.vector( piX ) ) nAlph <- 2
  else nAlph <- dim(piX)[[2]]+1
  Salph <- ( A - piMLE ) * cbind(1,piX)
  piXtpiX <- t( apply( cbind(1,piX), 1, function(x) outer(x,x) ) )
  VSalph <- matrix(colMeans( piMLE*(1-piMLE)*piXtpiX ), nAlph)
  ViSalph <- qr.solve( VSalph, t(Salph) )
  fitY <- function( eta ) {
    etaX <- c( cbind( 1, phiX ) %*% eta )
    if ( !YisBin ) return ( etaX )
    return( c( 1/( 1 + exp( -etaX ) ) ) )
  }
  if ( is.null(initEta) ) initEta <- coef(glm(Y~phiX, family=famY, subset=A==1))
  else {
    if ( length(initEta) != nEta ) stop("íinitEtaí should be a vector of ", nEta)
    if ( !is.numeric(initEta) ) stop("íinitEtaí should be numeric")
  }
  
  # Step 1. Function bounded returns  ^OR, Y ^ , and  ^OR
  bounded <- function( Pi ) {
    etaB <- glm(Y~phiX,family=famY,weights=1/Pi,subset=A==1,control=ctrl)$coef
    phiB <- fitY( etaB )
    list( etaB=etaB, phiB=phiB, betaB=mean( phiB ) )
  }
  OR <- bounded( piMLE )
  phiOR <- OR$phiB; betaOR <- OR$betaB
  
  # Step 2. Calculate  ~ minimizing estimated variance of 
  eff <- function( initEta, BETA ) {
    effObj <- function( eta, BETA ) {
      phihat = fitY( eta )
      psi <- YAdPi - BETA - ( AdPi - 1 )*phihat
      covar <- colMeans( psi * Salph )
      sum( ( psi - c( covar %*% ViSalph ) )^2 )
    }
    effgr <- function( eta, BETA ) {
      phihat = fitY( eta )
      psi <- YAdPi - BETA - ( AdPi - 1 )*phihat
      covar <- colMeans( psi * Salph )
      CF = 2 * ( psi - c( covar %*% ViSalph ) )
      EX = c( exp( - cbind(1,phiX) %*% eta ) )
      dEx <- 1 / ( 1/EX + 2 + EX )
      gr1 <- - ( AdPi - 1 ) * dEx * phiX1
      prj <- function(i) colMeans(-(AdPi-1)*dEx*phiX1[,i]*Salph)%*%ViSalph
      gr2 <- sapply( 1:nEta, prj )
      colSums( CF * ( gr1 + gr2 ) )
    }
    optim( initEta, effObj, gr=effgr, BETA=BETA, method="BFGS" )$par
  }
  effcont <- function( BETA ) {
    B <- ( AdPi - 1 ) * phiX1
    newY <- YAdPi - BETA
    newY <- newY - c( colMeans( newY*Salph ) %*% ViSalph )
    newX2 <- sapply( 1:nEta, function(i) colMeans(B[,i]*Salph)%*%ViSalph )
    newX <- B - newX2
    coef( glm( newY ~ 0 + newX ) )
  }
  if ( YisBin ) etaE <- eff( initEta, betaOR )
  else etaE <- effcont( betaOR )
  
  # Step 3. Fit extended model for missingness and repeat step 1 with Ötted missingness
  # probability from extended missingness model.
  phiE <- fitY( etaE )
  betaE = mean( YAdPi - ( AdPi - 1 ) * phiE )
  V1 = 1/piMLE * ( phiOR - betaOR )
  V2 = 1/piMLE * ( phiE - betaE )
  extX <- cbind( piX, V1, V2 )
  piext = glm( A ~ extX, family=binomial )$fitted
  DR <- bounded( piext )
  list( betaDR=DR$betaB, etaDR=DR$etaB )
}