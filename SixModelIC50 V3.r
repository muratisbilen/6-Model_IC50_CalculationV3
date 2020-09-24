tryCatch({
	library(car)
	}, error=function(e){
		install.packages("car",repos="http://cran.pau.edu.tr/")
		library(car)
	}
)

setwd("inputFiles")
inputFiles=try(shell("dir /B ",intern=T,wait=T))
setwd("../")

sigmoid3 = function(params, x) {
	params[1]+((params[2]-params[1])/(1+10^(x-params[3])))
}
sigmoid3B = function(params, x) {
	0+((params[1]-0)/(1+10^(x-params[2])))
}
sigmoid3T = function(params, x) {
	params[1]+((100-params[1])/(1+10^(x-params[2])))
}

sigmoid4 = function(params, x) {
	params[1]+((params[2]-params[1])/(1+10^((params[3]-x)*params[4])))
}
sigmoid4B = function(params, x) {
	0+((params[1]-0)/(1+10^((params[2]-x)*params[3])))
}
sigmoid4T = function(params, x) {
	params[1]+((100-params[1])/(1+10^((params[2]-x)*params[3])))
}

mat = c("Model","Sample","log(IC50)","Std. Err of Model","Std. Err of log(EC50)",
	"log(EC50)","Activity Area","Amax","log(IC90)","log(IC95)")

count = 0
pdf("Fitted Curves.pdf",width=14)
for(files in inputFiles){
	count = count+1
	setwd("inputFiles")
	d=read.delim(files, sep="\t",header=T)
	setwd("../")
	attach(d)
	
	l3="3-Parameter"
	l3B="3-Parameter Bottom 0"
	l3T="3-Parameter Top 100"
	l4="4-Parameter"
	l4B="4-Parameter Bottom 0"
	l4T="4-Parameter Top 100"

	x=log10(d[which(d[,1]!=0),1])
	y=d[which(d[,1]!=0),-1]
	if(is.vector(y)){
		y = cbind(y,y)
	}
	ym = apply(y,1,mean,na.rm=T)
	x2 = seq(min(x,na.rm=T),max(x,na.rm=T),0.01)
	
	l = c()
	for(i in 1:nrow(y)){
		if(length(y[i,which(!is.na(y[i,]))])>1){
			l = c(l,sd(y[i,],na.rm=T)/sqrt(length(which(!is.na(y[i,])))))
		}else{
			l = c(l,0)
		}
	}
	
	ylmax = max(ym+l)
	ylmin = min(ym-l)
	
	par(mfrow = c(2,3))
	tryCatch({
		fit3 = nls(ym~Ab+(At-Ab)/(1+10^(x-EC50)), start=list(Ab = min(ym), At=max(ym), EC50=0))
		params3 = coef(fit3)
		ic3=coef(fit3)[3]+log10((coef(fit3)[2]-coef(fit3)[1])/(50-coef(fit3)[1])-1)
		ic3.90=coef(fit3)[3]+log10((coef(fit3)[2]-coef(fit3)[1])/(10-coef(fit3)[1])-1)
		ic3.95=coef(fit3)[3]+log10((coef(fit3)[2]-coef(fit3)[1])/(5-coef(fit3)[1])-1)
		icse3 = summary(fit3)$sigma
		yp3 = sigmoid3(params3,x2)
		aa3 = sum(0.01*(params3[2]-yp3))
		amax3 = params3[2]-params3[1]
		ecse3 = summary(fit3)$coefficients[3,2]
		
		
		plot(x2,yp3,type="l",lwd=4,col="red",ylim=c(min(ylmin,yp3),max(ylmax,yp3)),
			xlab="Concentration (log µM)", ylab="Growth (%)",
			main=paste(substr(files,1,nchar(files)-4),"\n",l3,"\nError: ",format.pval(icse3),sep=""))
		points(x,ym,pch=16,cex=1.5)
		for(i in 1:nrow(y)){
			lines(c(x[i],x[i]),c(ym[i]-l[i],ym[i]+l[i]),lwd=2)
			lines(c(x[i]-0.05,x[i]+0.05),c(ym[i]-l[i],ym[i]-l[i]),lwd=2)
			lines(c(x[i]-0.05,x[i]+0.05),c(ym[i]+l[i],ym[i]+l[i]),lwd=2)
		}
	}, error=function(e){
		ic3<<-NA;ic3.90<<-NA;ic3.95<<-NA;icse3<<-NA;aa3<<-NA;amax3<<-NA;yp3<<-0;
		params3 <<- c(NA,NA,NA); ecse3 <<- NA
		plot(x,ym,pch=16,cex=1.5,ylim=c(min(ylmin,yp3),max(ylmax,yp3)),
			xlab="Concentration (log µM)", ylab="Growth (%)",
			main=paste(substr(files,1,nchar(files)-4),"\n",l3,sep=""))
		for(i in 1:nrow(y)){
			lines(c(x[i],x[i]),c(ym[i]-l[i],ym[i]+l[i]),lwd=2)
			lines(c(x[i]-0.05,x[i]+0.05),c(ym[i]-l[i],ym[i]-l[i]),lwd=2)
			lines(c(x[i]-0.05,x[i]+0.05),c(ym[i]+l[i],ym[i]+l[i]),lwd=2)
		}
	})
	
	tryCatch({
		fit3B = nls(ym~0+(At-0)/(1+10^(x-EC50)), start=list(At=max(ym), EC50=0))
		params3B = coef(fit3B)
		ic3B=coef(fit3B)[2]+log10((coef(fit3B)[1]-0)/(50-0)-1)
		ic3B.90=coef(fit3B)[2]+log10((coef(fit3B)[1]-0)/(10-0)-1)
		ic3B.95=coef(fit3B)[2]+log10((coef(fit3B)[1]-0)/(5-0)-1)
		icse3B = summary(fit3B)$sigma
		yp3B = sigmoid3B(params3B,x2)
		aa3B = sum(0.01*(params3B[1]-yp3B))
		amax3B = params3B[1]-0
		ecse3B = summary(fit3B)$coefficients[2,2]
		
		plot(x2,yp3B,type="l",lwd=4,col="red",ylim=c(min(ylmin,yp3B),max(ylmax,yp3B)),
			xlab="Concentration (log µM)", ylab="Growth (%)",
			main=paste(substr(files,1,nchar(files)-4),"\n",l3B,"\nError: ",format.pval(icse3B),sep=""))
		points(x,ym,pch=16,cex=1.5)
		for(i in 1:nrow(y)){
			lines(c(x[i],x[i]),c(ym[i]-l[i],ym[i]+l[i]),lwd=2)
			lines(c(x[i]-0.05,x[i]+0.05),c(ym[i]-l[i],ym[i]-l[i]),lwd=2)
			lines(c(x[i]-0.05,x[i]+0.05),c(ym[i]+l[i],ym[i]+l[i]),lwd=2)
		}
	}, error=function(e){
		ic3B<<-NA;ic3B.90<<-NA;ic3B.95<<-NA;icse3B<<-NA;aa3B<<-NA;amax3B<<-NA;yp3B<<-0;
		params3B <<- c(NA,NA);ecse3B <<- NA
		plot(x,ym,pch=16,cex=1.5,ylim=c(min(ylmin,yp3B),max(ylmax,yp3B)),
			xlab="Concentration (log µM)", ylab="Growth (%)",
			main=paste(substr(files,1,nchar(files)-4),"\n",l3B,sep=""))
		for(i in 1:nrow(y)){
			lines(c(x[i],x[i]),c(ym[i]-l[i],ym[i]+l[i]),lwd=2)
			lines(c(x[i]-0.05,x[i]+0.05),c(ym[i]-l[i],ym[i]-l[i]),lwd=2)
			lines(c(x[i]-0.05,x[i]+0.05),c(ym[i]+l[i],ym[i]+l[i]),lwd=2)
		}
	})

	tryCatch({
		fit3T = nls(ym~Ab+(100-Ab)/(1+10^(x-EC50)), start=list(Ab=min(ym), EC50=0))
		params3T = coef(fit3T)
		ic3T=coef(fit3T)[2]+log10((100-coef(fit3T)[1])/(50-coef(fit3T)[1])-1)
		ic3T.90=coef(fit3T)[2]+log10((100-coef(fit3T)[1])/(10-coef(fit3T)[1])-1)
		ic3T.95=coef(fit3T)[2]+log10((100-coef(fit3T)[1])/(5-coef(fit3T)[1])-1)
		icse3T = summary(fit3T)$sigma
		yp3T = sigmoid3T(params3T,x2)
		aa3T = sum(0.01*(100-yp3T))
		amax3T = 100-params3T[1]
		ecse3T = summary(fit3T)$coefficients[2,2]
		
		plot(x2,yp3T,type="l",lwd=4,col="red",ylim=c(min(ylmin,yp3T),max(ylmax,yp3T)),
			xlab="Concentration (log µM)", ylab="Growth (%)",
			main=paste(substr(files,1,nchar(files)-4),"\n",l3T,"\nError: ",format.pval(icse3T),sep=""))
		points(x,ym,pch=16,cex=1.5)
		for(i in 1:nrow(y)){
			lines(c(x[i],x[i]),c(ym[i]-l[i],ym[i]+l[i]),lwd=2)
			lines(c(x[i]-0.05,x[i]+0.05),c(ym[i]-l[i],ym[i]-l[i]),lwd=2)
			lines(c(x[i]-0.05,x[i]+0.05),c(ym[i]+l[i],ym[i]+l[i]),lwd=2)
		}
	}, error=function(e){
		ic3T<<-NA;ic3T.90<<-NA;ic3T.95<<-NA;icse3T<<-NA;aa3T<<-NA;amax3T<<-NA;yp3T<<-0;
		params3T <<- c(NA,NA);ecse3T <<- NA
		plot(x,ym,pch=16,cex=1.5,ylim=c(min(ylmin,yp3T),max(ylmax,yp3T)),
			xlab="Concentration (log µM)", ylab="Growth (%)",
			main=paste(substr(files,1,nchar(files)-4),"\n",l3T,sep=""))
		for(i in 1:nrow(y)){
			lines(c(x[i],x[i]),c(ym[i]-l[i],ym[i]+l[i]),lwd=2)
			lines(c(x[i]-0.05,x[i]+0.05),c(ym[i]-l[i],ym[i]-l[i]),lwd=2)
			lines(c(x[i]-0.05,x[i]+0.05),c(ym[i]+l[i],ym[i]+l[i]),lwd=2)
		}
	})
	tryCatch({
		fit4 = nls(ym~Ab+(At-Ab)/(1+10^((EC50-x)*H)), start=list(Ab = min(ym), At=max(ym), EC50=0,H=-1))
		params4 = coef(fit4)
		ic4=coef(fit4)[3]-log10((coef(fit4)[2]-coef(fit4)[1])/(50-coef(fit4)[1])-1)/coef(fit4)[4]
		ic4.90=coef(fit4)[3]-log10((coef(fit4)[2]-coef(fit4)[1])/(10-coef(fit4)[1])-1)/coef(fit4)[4]
		ic4.95=coef(fit4)[3]-log10((coef(fit4)[2]-coef(fit4)[1])/(5-coef(fit4)[1])-1)/coef(fit4)[4]
		icse4 = summary(fit4)$sigma
		yp4 = sigmoid4(params4,x2)
		aa4 = sum(0.01*(params4[2]-yp4))
		amax4 = params4[2]-params4[1]
		ecse4 = summary(fit4)$coefficients[3,2]
		
		plot(x2,yp4,type="l",lwd=4,col="red",ylim=c(min(ylmin,yp4),max(ylmax,yp4)),
			xlab="Concentration (log µM)", ylab="Growth (%)",
			main=paste(substr(files,1,nchar(files)-4),"\n",l4,"\nError: ",format.pval(icse4),sep=""))
		points(x,ym,pch=16,cex=1.5)
		for(i in 1:nrow(y)){
			lines(c(x[i],x[i]),c(ym[i]-l[i],ym[i]+l[i]),lwd=2)
			lines(c(x[i]-0.05,x[i]+0.05),c(ym[i]-l[i],ym[i]-l[i]),lwd=2)
			lines(c(x[i]-0.05,x[i]+0.05),c(ym[i]+l[i],ym[i]+l[i]),lwd=2)
		}
	}, error=function(e){
		ic4<<-NA;ic4.90<<-NA;ic4.95<<-NA;icse4<<-NA;aa4<<-NA;amax4<<-NA;yp4<<-0;
		params4 <<- c(NA,NA,NA,NA); ecse4 <<- NA
		plot(x,ym,pch=16,cex=1.5,ylim=c(min(ylmin,yp4),max(ylmax,yp4)),
			xlab="Concentration (log µM)", ylab="Growth (%)",
			main=paste(substr(files,1,nchar(files)-4),"\n",l4,sep=""))
		for(i in 1:nrow(y)){
			lines(c(x[i],x[i]),c(ym[i]-l[i],ym[i]+l[i]),lwd=2)
			lines(c(x[i]-0.05,x[i]+0.05),c(ym[i]-l[i],ym[i]-l[i]),lwd=2)
			lines(c(x[i]-0.05,x[i]+0.05),c(ym[i]+l[i],ym[i]+l[i]),lwd=2)
		}
	})
	
	tryCatch({
		fit4B = nls(ym~0+(At-0)/(1+10^((EC50-x)*H)), start=list(At=max(ym), EC50=0,H=-1))
		params4B = coef(fit4B)
		ic4B=coef(fit4B)[2]-log10((coef(fit4B)[1]-0)/(50-0)-1)/coef(fit4B)[3]
		ic4B.90=coef(fit4B)[2]-log10((coef(fit4B)[1]-0)/(10-0)-1)/coef(fit4B)[3]
		ic4B.95=coef(fit4B)[2]-log10((coef(fit4B)[1]-0)/(5-0)-1)/coef(fit4B)[3]
		icse4B = summary(fit4B)$sigma
		yp4B = sigmoid4B(params4B,x2)
		aa4B = sum(0.01*(params4B[1]-yp4B))
		amax4B = params4B[1]-0
		ecse4B = summary(fit4B)$coefficients[2,2]
		
		plot(x2,yp4B,type="l",lwd=4,col="red",ylim=c(min(ylmin,yp4B),max(ylmax,yp4B)),
			xlab="Concentration (log µM)", ylab="Growth (%)",
			main=paste(substr(files,1,nchar(files)-4),"\n",l4B,"\nError: ",format.pval(icse4B),sep=""))
		points(x,ym,pch=16,cex=1.5)
		for(i in 1:nrow(y)){
			lines(c(x[i],x[i]),c(ym[i]-l[i],ym[i]+l[i]),lwd=2)
			lines(c(x[i]-0.05,x[i]+0.05),c(ym[i]-l[i],ym[i]-l[i]),lwd=2)
			lines(c(x[i]-0.05,x[i]+0.05),c(ym[i]+l[i],ym[i]+l[i]),lwd=2)
		}
	}, error=function(e){
		ic4B<<-NA;ic4B.90<<-NA;ic4B.95<<-NA;icse4B<<-NA;aa4B<<-NA;amax4B<<-NA;yp4B<<-0;
		params4B <<- c(NA,NA,NA); ecse4B <<- NA
		plot(x,ym,pch=16,cex=1.5,ylim=c(min(ylmin,yp4B),max(ylmax,yp4B)),
			xlab="Concentration (log µM)", ylab="Growth (%)",
			main=paste(substr(files,1,nchar(files)-4),"\n",l4B,sep=""))
		for(i in 1:nrow(y)){
			lines(c(x[i],x[i]),c(ym[i]-l[i],ym[i]+l[i]),lwd=2)
			lines(c(x[i]-0.05,x[i]+0.05),c(ym[i]-l[i],ym[i]-l[i]),lwd=2)
			lines(c(x[i]-0.05,x[i]+0.05),c(ym[i]+l[i],ym[i]+l[i]),lwd=2)
		}
	})
	
	tryCatch({	
		fit4T = nls(ym~Ab+(100-Ab)/(1+10^((EC50-x)*H)), start=list(Ab = min(ym), EC50=0,H=-1))
		params4T = coef(fit4T)
		ic4T=coef(fit4T)[2]-log10((100-coef(fit4T)[1])/(50-coef(fit4T)[1])-1)/coef(fit4T)[3]
		ic4T.90=coef(fit4T)[2]-log10((100-coef(fit4T)[1])/(10-coef(fit4T)[1])-1)/coef(fit4T)[3]
		ic4T.95=coef(fit4T)[2]-log10((100-coef(fit4T)[1])/(5-coef(fit4T)[1])-1)/coef(fit4T)[3]
		icse4T = summary(fit4T)$sigma
		yp4T = sigmoid4T(params4T,x2)
		aa4T = sum(0.01*(100-yp4T))
		amax4T = 100-params4T[1]
		ecse4T = summary(fit4T)$coefficients[2,2]
		
		plot(x2,yp4T,type="l",lwd=4,col="red",ylim=c(min(ylmin,yp4T),max(ylmax,yp4T)),
			xlab="Concentration (log µM)", ylab="Growth (%)",
			main=paste(substr(files,1,nchar(files)-4),"\n",l4T,"\nError: ",format.pval(icse4T),sep=""))
		points(x,ym,pch=16,cex=1.5)
		for(i in 1:nrow(y)){
			lines(c(x[i],x[i]),c(ym[i]-l[i],ym[i]+l[i]),lwd=2)
			lines(c(x[i]-0.05,x[i]+0.05),c(ym[i]-l[i],ym[i]-l[i]),lwd=2)
			lines(c(x[i]-0.05,x[i]+0.05),c(ym[i]+l[i],ym[i]+l[i]),lwd=2)
		}
	}, error=function(e){
		ic4T<<-NA;ic4T.90<<-NA;ic4T.95<<-NA;icse4T<<-NA;aa4T<<-NA;amax4T<<-NA;yp4T<<-0;
		params4T <<- c(NA,NA,NA); ecse4T <<- NA
		plot(x,ym,pch=16,cex=1.5,ylim=c(min(ylmin,yp4T),max(ylmax,yp4T)),
			xlab="Concentration (log µM)", ylab="Growth (%)",
			main=paste(substr(files,1,nchar(files)-4),"\n",l4T,sep=""))
		for(i in 1:nrow(y)){
			lines(c(x[i],x[i]),c(ym[i]-l[i],ym[i]+l[i]),lwd=2)
			lines(c(x[i]-0.05,x[i]+0.05),c(ym[i]-l[i],ym[i]-l[i]),lwd=2)
			lines(c(x[i]-0.05,x[i]+0.05),c(ym[i]+l[i],ym[i]+l[i]),lwd=2)
		}
	})
	
	r_3 = c(l3,substr(files,1,nchar(files)-4),ic3,icse3,ecse3,params3[3],aa3,amax3,
		ic3.90,ic3.95)
	r_3B = c(l3B,substr(files,1,nchar(files)-4),ic3B,icse3B,ecse3B,params3B[2],aa3B,
		amax3B,ic3B.90,ic3B.95)
	r_3T = c(l3T,substr(files,1,nchar(files)-4),ic3T,icse3T,ecse3T,params3T[2],aa3T,
		amax3T,ic3T.90,ic3T.95)
	r_4 = c(l4,substr(files,1,nchar(files)-4),ic4,icse4,ecse4,params4[3],aa4,amax4,
		ic4.90,ic4.95)
	r_4B = c(l4B,substr(files,1,nchar(files)-4),ic4B,icse4B,ecse4B,params4B[2],aa4B,
		amax4B,ic4B.90,ic4B.95)
	r_4T = c(l4T,substr(files,1,nchar(files)-4),ic4T,icse4T,ecse4T,params4T[2],aa4T,
		amax4T,ic4T.90,ic4T.95)
	
	
	mat = rbind(mat,r_3,r_3B,r_3T,r_4,r_4B,r_4T)
	detach(d)
	print(paste(count,length(inputFiles),sep="/"))
	flush.console()
}

dev.off()
write(file="DrugData.txt",t(mat),ncol=ncol(mat),sep="\t")
