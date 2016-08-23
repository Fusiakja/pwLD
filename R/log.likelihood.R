`log.likelihood` <-
function(case, code)
{
	
	 if(case == paste(code[1], "-", code[1], sep="")) return("2*log(Poo)")							# 0 - 0
	 if(case == paste(code[1], "-", code[2], sep="")) return("(log(2)+log(Poo)+log(po.-Poo))")				# 0 - 1	
         if(case == paste(code[1], "-", code[3], sep="")) return("2*log(po.-Poo)")					   	# 0 - 2
	 if(case == paste(code[2], "-", code[1], sep="")) return("(log(2)+log(Poo)+log(p.o-Poo))")				# 1 - 0
	 if(case == paste(code[2], "-", code[2], sep="")) return("log(2*Poo*(1-po.-p.o+Poo) + 2*(po.-Poo)*(p.o-Poo))")      	# 1 - 1
	 if(case == paste(code[2], "-", code[3], sep="")) return("(log(2)+log(po.-Poo)+log(1-po.-p.o+Poo))")			# 1 - 2
	 if(case == paste(code[3], "-", code[1], sep="")) return("(2*log(p.o-Poo))")				        	# 2 - 0
	 if(case == paste(code[3], "-", code[2], sep="")) return("(log(2)+log(p.o-Poo)+log(1-po.-p.o+Poo))")			# 2 - 1
	 if(case == paste(code[3], "-", code[3], sep="")) return("2*log(1-po.-p.o+Poo)")						# 2 - 2
       	 if(case == paste(code[1], "-", code[4], sep="")) return("2*log(po.)")							# 0 - -1
       	 if(case == paste(code[2], "-", code[4], sep="")) return("(log(2)+log(po.)+log(1-po.))")			              	# 1 - -1
	 if(case == paste(code[3], "-", code[4], sep="")) return("2*log(1-po.)")							# 2 - -1
	 if(case == paste(code[4], "-", code[1], sep="")) return("2*log(p.o)") 							# -1 - 0 
  	 if(case == paste(code[4], "-", code[2], sep="")) return("(log(2)+log(p.o)+log(1-p.o))")			         	# -1 - 1
	 if(case == paste(code[4], "-", code[3], sep="")) return("2*log(1-p.o)")							# -1 - 2
	 if(case == paste(code[4], "-", code[4], sep="")) return("0")				                        	# -1 - -1

}

