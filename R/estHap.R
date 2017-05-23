`estHap` <-
function( genoFreq, tol, digits=12, Neff)
{
	res=.C("estimateHaploFreq", genoFreq=as.double(t(genoFreq)), haploFreq=as.double(matrix(0, nrow=3, ncol=3)), tol=as.double(tol), digits=as.integer(digits), pooHat=as.double(rep(0,3)), PACKAGE="pwLD" )
  # genoCount <- matrix(rep(0,12),nrow = 4, ncol = 4)
  # genoCount[1,1] <- Neff*res$haploFreq[1]^2
  # genoCount[1,2] <- 2*Neff*res$haploFreq[1]*res$haploFreq[2]
  # genoCount[1,3] <- Neff*res$haploFreq[2]^2
  # genoCount[2,1] <- 2*Neff*res$haploFreq[1]*res$haploFreq[4]
  # genoCount[2,2] <- 2*Neff*res$haploFreq[2]*res$haploFreq[4]+2*Neff*res$haploFreq[1]*res$haploFreq[5]
  # genoCount[2,3] <- 2*Neff*res$haploFreq[2]*res$haploFreq[5]
  # genoCount[3,1] <- Neff*res$haploFreq[4]^2
  # genoCount[3,2] <- 2*Neff*res$haploFreq[4]*res$haploFreq[5]
  # genoCount[3,3] <- Neff*res$haploFreq[5]^2
 #  if (genoFreq[2,2]==0) {
 # 
 #    haploFreq <- matrix(rep(0,9), nrow = 3, ncol = 3)
 #    haploFreq[1,1] <- 2*genoFreq[1,1]+genoFreq[1,2]+genoFreq[2,1]
 #    haploFreq[2,1] <- 2*genoFreq[1,3]+genoFreq[1,2]+genoFreq[2,3]
 #    haploFreq[1,2] <- 2*genoFreq[3,1]+genoFreq[3,2]+genoFreq[2,1]
 #    haploFreq[2,2] <- 2*genoFreq[3,3]+genoFreq[3,2]+genoFreq[2,3]
 # 
 #  }
 #  else
 #  {
 #    Neff <- sum(genoFreq)
 #    haploFreq <- matrix(rep(0,9), nrow = 3, ncol = 3)
 #    haploFreq[1,1] <- 2*genoFreq[1,1]+genoFreq[1,2]+genoFreq[2,1]
 #    haploFreq[2,1] <- 2*genoFreq[1,3]+genoFreq[1,2]+genoFreq[2,3]
 #    haploFreq[1,2] <- 2*genoFreq[3,1]+genoFreq[3,2]+genoFreq[2,1]
 #    haploFreq[2,2] <- 2*genoFreq[3,3]+genoFreq[3,2]+genoFreq[2,3]
 #    p <- (2*(genoFreq[1,1]+genoFreq[1,2]+genoFreq[1,3])+genoFreq[2,1]+genoFreq[2,2]+genoFreq[2,3])/(2*Neff)
 #    q <- (2*(genoFreq[1,1]+genoFreq[2,1]+genoFreq[3,1])+genoFreq[1,2]+genoFreq[2,2]+genoFreq[3,2])/(2*Neff)
 # 
 #   hwe <- p*p+2*p*q+q*q
 #    a <- 4*Neff
 #    b <- 2*Neff*(1-2*p-2*q)-2*(2*haploFreq[1,1]+haploFreq[1,2]+haploFreq[2,1])-haploFreq[2,2]
 #    c <- 2*Neff*p*q-(2*haploFreq[1,1]+haploFreq[1,2]+haploFreq[2,1])*(1-2*p-2*q)-haploFreq[2,2]*(1-p-q)
 #    d <- -(2*haploFreq[1,1]+haploFreq[1,2]+haploFreq[2,1])*p*q
 # 
 #    x_n <- -(b/(3*a))
 #    d_square <- (b*b-3*a*c)/(9*a*a)
 #    h_square <- (4*a*a*d_square*d_square*d_square)
 #    y_n <- a*x_n*x_n*x_n+b*x_n*x_n+c*x_n+d
 # 
 #    delta <- ((y_n*y_n)-h_square)
 # 
 #    if ((y_n*y_n) > h_square)
 #    {
 # 
 #      alphas <- x_n+((1/(2*a))*(-y_n+sqrt(delta)))^(1/3)+((1/(2*a))*(-y_n-sqrt(delta)))^(1/3)
 #      betas <- -1
 #      gammas <- -1
 #    }
 #    else if ((y_n*y_n) == h_square)
 #    {
 #      mu <- (y_n/(2*a))^(1/3)
 # 
 #      alphas <- x_n+mu
 #      betas <- x_n+mu
 #      gammas <- x_n-2*mu
 #    }
 #    else if ((y_n*y_n) < h_square)
 #    {
 # 
 #      theta <- acos(-y_n/sqrt(h_square))/3
 # 
 #      alphas <- x_n+2*sqrt(d_square)*cos(theta)
 #      betas <- x_n+2*sqrt(d_square)*cos((2*pi/3)+theta)
 #      gammas <- x_n+2*sqrt(d_square)*cos((4*pi/3)+theta)
 #    }
 # 
 # }

	output <- list()
	output[[1]] <- matrix(res$haploFreq, nrow=3, byrow=T, dimnames=list(c("0", "1", "Sum"), c("0", "1", "Sum")))
	output[[2]] <- res$pooHat

	names(output) <- c("table", "solution")
	return( output )
}

