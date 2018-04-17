
dnorm_product <- function(mu1,sig1,mu2,sig2){
	mu <- (mu1*sig2+mu2*sig1)/(sig2+sig1)
	sig <- (sig1*sig2)/(sig1+sig2)
	const <- (exp(-((mu1-mu2)^2)/(2*(sig1+sig2))))/sqrt(2*pi*(sig1+sig2))
	return(c(mu,sig,const))
}


##We supposed that mu1 is always inferior to mu2
overlap_n1_n2<- function(mu1,sig1,mu2,sig2,fac_sig=2){


	bmin <- mu2-fac_sig*sqrt(sig2)
	bmax <- mu1+fac_sig*sqrt(sig1)
		##We check that any of the masses is superior to this hsit.

	if(bmax < bmin){
		# cat("res",0,"\n")
		return(0)
	}

	area1 <- pnorm(bmax,mu1,sqrt(sig1))-pnorm(bmin,mu1,sqrt(sig1))
	area2 <- pnorm(bmax,mu2,sqrt(sig2))-pnorm(bmin,mu2,sqrt(sig2))
	mval <- max(area1,area2)
	return(max(area1,area2))

	# npar <- dnorm_product(mu1,sig1,mu2,sig2)
	# area <- (pnorm(bmax,npar[1],sqrt(npar[2]))-pnorm(bmin,npar[1],sqrt(npar[2])))*npar[3]
	# return((pnorm(bmax,npar[1],sqrt(npar[2]))-pnorm(bmin,npar[1],sqrt(npar[2])))*npar[3])
}

mergeGaussian <- function(mu,sig,fac_sig=2){
	inter <- fac_sig*sqrt(sig)
	bmin <- mu-inter
	bmax <- mu+inter
	inter <- c(min(bmin),max(bmax))
	sd_v <- diff(inter)/(2*fac_sig)
	mean_v <- sum(inter)/2
	return(c(mean_v,sd_v^2))
}



###A naive algoirthm which split the gaussian according to the dataset.
gaussianMerging<- function(mu,sig,alpha=0.01,fac_sig=3){
	mu <- c(mu,1000)
	sig <- c(sig,0.01)
	start <- NA
	stop <- NULL
	res <- list()
	res$mu <- numeric(length(mu))
	res$sig <- numeric(length(mu))
	res$fused <- logical(length(mu))
	nres <- 0
	nmerge <- 0
	if(length(mu)==1)
		return(mu_si)
	for(i in 1:(length(sig)-1)){
		###We get the value of the parameter
		mu1 <- mu[i];sig1 <- sig[i]
		mu2 <- mu[i+1];sig2 <- sig[i+1]
		ov <- overlap_n1_n2(mu1,sig1,mu2,sig2,fac_sig = fac_sig)

		if(ov>alpha & is.na(start)){
			start <- i
		}else if(ov<alpha & !is.na(start)){
			stop <- i
			ng <- mergeGaussian(mu[start:stop],sig[start:stop],fac_sig)
			nres <- nres+1
			nmerge <- nmerge+1
			res$mu[nres]<-ng[1]
			res$sig[nres]<-ng[2]
			res$fused[nres]<-TRUE

			start <- NA
		}else if(ov < alpha & is.na(start)){
			nres <- nres+1
			res$mu[nres]<-mu[i]
			res$sig[nres]<-sig[i]
			res$fused[nres]<-FALSE
		}
	}
	cat("Passed from",length(mu),"to",nres,"points.\n")
	res$mu <- res$mu[1:nres]
	res$sig <- res$sig[1:nres]
	res$fused<-res$fused[1:nres]
	return(res)
}
