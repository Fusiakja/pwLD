`cardano.complex` <-
function( coeff )
{
	if(coeff[1] == 0) stop("not a cubic equation!!")

	output <- vector("list", 2)
	names(output) <- c("D", "x")

	a <- coeff[1]; b <- coeff[2]; c <- coeff[3]; d <- coeff[4]


	# some special cases
        # linear function
	if(c==0 & d==0)
	{
		x <- -b/a

		output[1] <- NA
		output[[2]] <- c(0,0,x) 	
		
		return(output)
	}
	

        # quadratic function
	if(d == 0)
	{
		x1 <- 0

		dis <- (b/a/2)^2-(c/a)
		if(dis == 0)
		{
			output[1] <- NA
			output[2] <- x1

			return(output)
		}
		if(dis > 0)	
		{	 
			x2 = -(b/a/2) + sqrt(dis)
		 	x3 = -(b/a/2) - sqrt(dis)

			output[1] <- NA
			output[[2]] <- c(x1,x2,x3)

			return(output)
		}
				
	}

	# reduced form of a cubic equation: x^3 + px +q = 0
	p <- ((3*a*c - b^2) / (3*a^2))/3
	q <- ((2*b^3)/(27*a^3) - (b*c)/(3*a^2) + (d/a))/2

	D <- p^3 + q^2

	# calculation of roots within the set of complex numbers
	u1 <- as.complex(-q+ sqrt(as.complex(D)))
	if(Im(u1)==0 & Re(u1) < 0) u <- -abs(u1)^(1/3)
	else u <- u1^(1/3)
 
	v1 <- as.complex(-q - sqrt(as.complex(D)))
	if(Im(v1)==0 & Re(v1) < 0) v <- -abs(v1)^(1/3)
	else v <- v1^(1/3)

	#####
	f1 <- 0.5*complex(1,-1,sqrt(3))
	f2 <- 0.5*complex(1,-1,-sqrt(3))
	#####

	y <- c(u+v, f1*u + f2*v, f2*u + f1*v)
	

	x <- y - b/(3*a)

	output[1] <- D
	output[[2]] <- x

return(output)

}

