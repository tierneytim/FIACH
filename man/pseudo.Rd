\name{pseudo}
\alias{pseudo}
\title{
Pseudoinverse.
}
\description{
This function computes the Pseudoinverse of a design matrix. If a matrix of dependent variables is also supplied the betas will also be returned. Optionally the residuals may be returned as well.}
\usage{
pseudo(x,y=NULL,residuals=FALSE,keepMean=FALSE,includeIntercept=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
x must be a numeric matrix.
}
  \item{y}{
y must be a numeric matrix.
}
 \item{residuals}{
residuals must be either TRUE or FALSE.
}

 \item{keepMean}{
keepMean must be either TRUE or FALSE. This argument specifies whether residuals should be zero centered or retain their mean. 
}
 \item{includeIntercept}{
includeIntercept must be either TRUE or FALSE. This argument specifies whether a column of ones should be included in x. 
}
}
\value{
Returns a Matrix containing either the pseudoinverse, linear regression coefficients or residuals.
}

\examples{

data(trees, package="datasets")

## with intercept
pinv<-pseudo(x=log(trees$Girth),includeIntercept=TRUE)
## without intercept
pinv<-pseudo(x=log(trees$Girth),includeIntercept=FALSE)
##coefficients
coef<-pseudo(x=log(trees$Girth),y=log(trees$Volume))
## residuals
res<-pseudo(x=log(trees$Girth),y=log(trees$Volume),residuals=TRUE)

## standard model
mod<-lm(log(trees$Volume)~log(trees$Girth))
coef2<-coefficients(mod)
res2<-residuals(mod)
## equality

cbind(coef,coef2) ## same but order reversed
cbind(res,res2)   ## same


}
