#!/share/apps/R-3.6.1/bin/Rscript

# Provide TestName as first argument on the command line (or set it below)
# An input file should have name TestName.input.txt and tab-separated columns:
# Dist	NSims	NCont	NCase	ContMean	CaseMean	ContVar	CaseVar
# Each row defines one set of simulations. 
# The Dist variable must be either Normal or Poisson
# If Poisson is chosen the ContVar and CaseVar variables are ignored
# See example input files provided in the git repository

# David Curtis 2023 d.curtis@ucl.ac.uk

library(ggplot2)

wd="/Users/dave_000/OneDrive/sharedseq/TestTTest"
wd="/Users/Dave/OneDrive/sharedseq/TestTTest"
# setwd(wd)

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=1) {
  print("Provide TestName argument")
  quit()
}

TestName=args[1]
InputFile=sprintf("%s.input.txt",TestName)
ResultsFile=sprintf("%s.output.txt",TestName)

Input=data.frame(read.table(InputFile,header=TRUE,sep="\t",stringsAsFactors=FALSE))

qqSLP<-function(SLP,bottom,top,Filename) {
	SLP=SLP[order(-SLP)]
	rankSLP <- rank(SLP, na.last=TRUE, ties.method="first")
	nelem=length(SLP)
	midrank <- (nelem+1)/2
	rMinusMr <- rankSLP-midrank
	absDiff <- abs(rMinusMr/midrank)
	pVal <- 1-absDiff
	logP <- log10(pVal)
	eSLP <- sign(rankSLP - midrank) * -logP
	ppi=600
	png(Filename,width=6*ppi, height=6*ppi, res=ppi)

	toPlot=data.frame(matrix(ncol=2,nrow=length(SLP)))
	toPlot[,1]=SLP
	toPlot[,2]=eSLP
	colnames(toPlot)=c("SLP","eSLP")
	myplot=ggplot(toPlot,aes_q(x=as.name("eSLP"),y=as.name("SLP")))+geom_point(size=1)+ theme_bw() + 
		geom_hline(yintercept=0,size=1.0) +
		geom_vline(xintercept=0,size=1.0) +
		theme(panel.grid.major=element_line(colour = "black",size=0.25)) +
		scale_x_continuous(breaks = seq(2*floor(bottom/2),2*ceiling(top/2),by =2),minor_breaks=NULL,limits=c(2*floor(bottom/2),2*ceiling(top/2))) +
		scale_y_continuous(breaks = seq(2*floor(bottom/2),2*ceiling(top/2),by =2),minor_breaks=NULL,limits=c(2*floor(bottom/2),2*ceiling(top/2))) 
	print(myplot)
	dev.off()
}

# Function to calculate the likelihood ratio chi-squared for a comparison of two means, using logistic regression
LR.chisq <- function(x, y) {
		Data=data.frame(matrix(ncol=2,nrow=length(x)+length(y)))
		colnames(Data)=c("Outcome","Var")
		Data$Outcome[1:length(x)]=0
		Data$Var[1:length(x)]=x
		Data$Outcome[(length(x)+1):(length(x)+length(y))]=1
		Data$Var[(length(x)+1):(length(x)+length(y))]=y
		if (sum(Data$Var)==0) {
			ch2=0
		} else {
		m=glm(as.formula("Outcome~ 1"),data=Data,family="binomial")
		LL0=as.numeric(logLik(m))
		m=glm(as.formula("Outcome~ Var"),data=Data,family="binomial")
		LL1=as.numeric(logLik(m))
		ch2=2*(LL1-LL0)
		}
  return(ch2)
}

Results=data.frame(matrix(ncol=3,nrow=nrow(Input)))
colnames(Results)=c("UneqVarTT","EqVarTT","LR")
MeanResults=Results
colnames(MeanResults)=sprintf("%sMean",colnames(Results))
UnderResults05=Results
colnames(UnderResults05)=sprintf("%sUnder05",colnames(Results))
OverResults05=Results
colnames(OverResults05)=sprintf("%sOver05",colnames(Results))
UnderResults4=Results
colnames(UnderResults4)=sprintf("%sUnder00001",colnames(Results))
OverResults4=Results
colnames(OverResults4)=sprintf("%sOver00001",colnames(Results))

for (r in 1:nrow(Input)) {
	SimResults=data.frame(matrix(ncol=ncol(Results),nrow=Input$NSims[r]))
	colnames(SimResults)=colnames(Results)
	for (s in 1:Input$NSims[r]) {
		Seed=s+(r-1)*Input$NSims
		set.seed(Seed)
		if (Input$Dist[r]=="Poisson") {
			Cases=rpois(Input$NCase[r],Input$CaseMean[r])
			Controls=rpois(Input$NCont[r],Input$ContMean[r])
		} else if (Input$Dist[r]=="Normal") {
			Cases=rnorm(Input$NCase[r],Input$CaseMean[r],sqrt(Input$CaseVar[r]))
			Controls=rnorm(Input$NCont[r],Input$ContMean[r],sqrt(Input$ContVar[r]))
		} else {
			print(sprintf("Unrecognised distribution: ",Input$Dist[r]))
			quit()
		}
		tt=t.test(Controls,Cases,var.equal=FALSE)
		p=tt$p.value
		SimResults$UneqVarTT[s]=log10(p)*as.numeric(sign(tt$estimate[1]-tt$estimate[2]))
		tt=t.test(Controls,Cases,var.equal=TRUE)
		p=tt$p.value
		SimResults$EqVarTT[s]=log10(p)*as.numeric(sign(tt$estimate[1]-tt$estimate[2]))
		ch2=LR.chisq(Controls,Cases)
		p=1-pchisq(q=ch2, df=1)
		SimResults$LR[s]=log10(p)*sign(mean(Controls)-mean(Cases))
	}
	write.table(SimResults,sprintf("%s.SimResults.%d.txt",TestName,r),row.names=FALSE,quote=FALSE,sep="\t")
	for (c in 1:ncol(SimResults)) {
		qqSLP(SimResults[,c],min(SimResults),max(SimResults),sprintf("%s.SimResults.%d.%s.png",TestName,r,colnames(SimResults)[c]))
	}
	for (c in 1:ncol(SimResults)) {
		UnderResults05[r,c]=length(which(SimResults[,c]<log10(0.05)))
		OverResults05[r,c]=length(which(SimResults[,c]>-log10(0.05)))
		UnderResults4[r,c]=length(which(SimResults[,c]< (-4)))
		OverResults4[r,c]=length(which(SimResults[,c]>4))
		MeanResults[r,c]=mean(SimResults[,c])
	}
}

Results=cbind(MeanResults,UnderResults05)
Results=cbind(Results,OverResults05)
Results=cbind(Results,UnderResults4)
Results=cbind(Results,OverResults4)

write.table(Results,ResultsFile,row.names=FALSE,quote=FALSE,sep="\t")


