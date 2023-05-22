library(Rcpp)

# solves t(beta)%*%xx%*%beta/2 - t(xy)%*%beta + fac||beta||_1 
#or ||X%*%bet-y||_2^2/2 + fac||beta||_1 with xx=crossprod(X), xy=crossprod(x,y)
# subject to t(cmat)%*%beta=0
cppFunction('NumericVector ConstrLassoC0(NumericVector betstart, NumericMatrix xx, NumericVector xy, NumericMatrix cmat, double fac, int maxiter, double tol){
int p = xx.nrow();
int m = cmat.ncol();
NumericVector bet(p);
NumericVector bet_old(p);
NumericVector ksi(m);
NumericVector ksi_old(m);
NumericMatrix xxcc(p,p);
double tmp;
double tmp2;
double fac2;
double myabs;
double mymax;
int mysgn;
int iter;
int iter2;
int k;
double eps;
double eps2;
LogicalVector nonzero(p);
IntegerVector ind(p);


for (int i=0; i<p; i++){
    bet[i] = betstart[i];
    nonzero[i] = xx(i,i)>0;
    ind[i] = i;
}

for (int i=0; i<m; i++){
    ksi[i] = 0;
}

for (int i=0; i<p; i++){
    for (int j=0; j<p; j++){
        xxcc(i,j) = xx(i,j);
        for (int ij=0; ij<m; ij++){
            xxcc(i,j) += cmat(i,ij)*cmat(j,ij);
        }
    }
}

iter2 = 1;
eps2 = 1;
while (eps2>tol & iter2<maxiter){
    iter = 1;
    eps = 1;
    while (eps>tol & iter<maxiter){
        for (int i=0; i<p; i++){
            bet_old[i] = bet[i];
        }
        std::random_shuffle(ind.begin(), ind.end());
        for (int i=0; i<p; i++){
            k = ind[i];
            if(nonzero[k]){
                fac2 = fac/xxcc(k,k);
                tmp = 0;
                for (int j=0; j<p; j++){
                    tmp += xxcc(k,j)*bet[j];
                }
                tmp = tmp - xxcc(k,k)*bet[k];
                tmp2 = 0;
                for (int j=0; j<m; j++){
                    tmp2 += cmat(k,j)*ksi[j];
                }
                tmp = (xy[k] - tmp2 - tmp)/xxcc(k,k);
                if (tmp>0){
                    mysgn = 1;
                    myabs = tmp;
                }else if (tmp<0){
                    mysgn = -1;
                    myabs = -tmp;
                }else{
                    mysgn = 0;
                    myabs = 0;
                }
                if (myabs > fac2){
                    bet[k] = mysgn*(myabs - fac2);
                }else{
                    bet[k] = 0;
                }
            }else{
                 bet[k] = 0;
            }
        }
        eps = 0;
        for (int i=0; i<p; i++){
            eps += pow(bet[i] - bet_old[i], 2.0);
        }
        eps = sqrt(eps);
        iter += 1;
    }
    for (int i=0; i<m; i++){
        ksi_old[i] = ksi[i];
        for (int j=0; j<p; j++){
            ksi[i] += cmat(j,i)*bet[j];
        }
    }
    eps2 = 0;
    for (int j=0; j<m; j++){
        eps2 += pow(ksi[j] - ksi_old[j], 2.0);
    }
    eps2 = sqrt(eps2);
    iter2 += 1;  
}
if(iter2==maxiter){
    for (int i=0; i<p; i++){
        bet[i] = 0;
    }
}
return bet;
}')




ConstrLasso0 <- function(y, x, C, lambda, betstart=NULL, maxiter=1000, tol=1e-8){
# solves ||y-X%*%beta||_2^2/(2n) + lam||beta||_1 subject to t(C)%*%beta=0
    C <- as.matrix(C)
    XX <- crossprod(x)
    Xy <- crossprod(x,y)
    fac <- length(y)*lambda
    if (is.null(betstart)) betstart <- rep(0, ncol(x))
    set.seed(0)
    bet <- ConstrLassoC0(betstart, XX, Xy, C, fac, maxiter, tol) 
    return(bet)
}

ConstrLasso <- function(y, x, C=NULL, lambda=NULL, nlam=20, intercept=TRUE, scaling=TRUE, maxiter=1000, tol=1e-8){
    n <- nrow(x)
    p <- ncol(x)
    if(is.null(C)) C <- matrix(0,p,1)
    
    if (intercept){
        y.mean <- mean(y)
        y <- y - y.mean
        x.mean <- colMeans(x)
        x <- scale(x, center=x.mean, scale=F) # centering before scaling
    }
    if (scaling){
        x.sd <- apply(x,2,sd)
        x <- scale(x, center=F, scale=x.sd)
        C <- C/x.sd
    }
    
    # get sequence of tuning parameter lambda
    if (is.null(lambda)){
        maxlam <- 2*max(abs(crossprod(x,y)/n))
        lambda <- exp(seq(from=log(maxlam), to=log(1e-4), length.out=nlam))
    }
    nlam <- length(lambda)

    C <- as.matrix(C)
    xx <- crossprod(x)
    xy <- crossprod(x,y)
    fac <- length(y)*lambda
    betstart <- rep(0, ncol(x))
    set.seed(0)
    bet <- mapply(ConstrLassoC0, fac=fac, MoreArgs=list(xy=xy, xx=xx, cmat=C, betstart=betstart, maxiter=maxiter, tol=tol))
    
    if (scaling){
        bet <- bet/x.sd
    }
    int <- rep(0, nlam)
    if (intercept){
        int <- y.mean - as.vector(x.mean%*%bet)
    }
    return(list(int=int, bet=bet, lambda=lambda))
}



ConstrLassoCrossVal <- function(y, x, C=NULL, lambda=NULL, nlam=20, intercept=TRUE, scaling=TRUE, nfolds=5, maxiter=1000, tol=1e-8){
    method <- "ConstrLasso"
    n <- nrow(x)
    p <- ncol(x)
    if(is.null(C)) C <- matrix(0,p,1)
    
    y.ori <- y
    x.ori <- x
    if (intercept){
        y.mean <- mean(y)
        y <- y - y.mean
        x.mean <- colMeans(x)
        x <- scale(x, center=x.mean, scale=F) # centering before scaling
    }
    if (scaling){
        x.sd <- apply(x,2,sd)
        x <- scale(x, center=F, scale=x.sd)
        C <- C/x.sd
    }
    
    # get sequence of tuning parameter lambda
    if (is.null(lambda)){
        maxlam <- 2*max(abs(crossprod(x,y)/n))
        lambda <- exp(seq(from=log(maxlam), to=log(1e-4), length.out=nlam))
    }
    nlam <- length(lambda)


    # run the cv folds
    err <- matrix(0, nlam, nfolds)
    set.seed(0)
    rnum <- sample.int(n)%%nfolds+1
    for (j in 1:nfolds){
        temp <- rnum==j
        print(j)
        ytrain <- y[!temp]
        xtrain <- x[!temp,]
        ytest <- y[temp]
        xtest <- x[temp,]
        res <- do.call(method, list(y=ytrain, x=xtrain, C=C, lambda=lambda, intercept=FALSE, scaling=FALSE, maxiter=maxiter, tol=tol))
        # with annotation
        err[,j] <- colMeans((xtest%*%res$bet-res$int-as.vector(ytest))^2)
    }
    cvm <- rowMeans(err)
    cvsd <- apply(err,1,sd)/sqrt(nfolds)
    # fit with all lambda
    res.fit <- do.call(method, list(y=y, x=x, C=C, lambda=lambda, intercept=FALSE, scaling=FALSE, maxiter=maxiter, tol=tol))
    bet <- res.fit$bet
    rownames(bet) <- colnames(xmat)
    if (scaling){
        bet <- bet/x.sd
    }
    int <- rep(0, nlam)
    if (intercept){
        int <- y.mean - as.vector(x.mean%*%bet)
    }  
    sel <- min(which(cvm <= min(cvm+cvsd)))
    bet.sel <- bet[,sel]
    int.sel <- int[sel]
    Rsq <- 1 - sum((int.sel+x.ori%*%bet.sel - y.ori)^2) / sum((y.ori-mean(y.ori))^2)
    return(list(lambda=lambda, int=int, bet=bet, cvm=cvm, cvsd=cvsd, 
        int.sel=int.sel, bet.sel=bet.sel, sel=sel, Rsq.sel=Rsq, 
        x=x.ori, y=y.ori))
}

