#!/share/apps/R-3.6.1/bin/Rscript

wd="/Users/dave_000/OneDrive/sharedseq/ttest"
wd="/Users/Dave/OneDrive/sharedseq/ttest"
setwd(wd)

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=1) {
  print("Provide TestNameargument")
  quit()
}

TestName=args[1]
InputFile=sprintf("%s.input.txt",TestName)
ResultsFile=sprintf("%s.output.txt",TestName)

Input=data.frame(read.table(InputFile,header=TRUE,sep="\t",stringsAsFactors=FALSE))

t.test.unequal.var <- function(x, y) {
  # Calculate the means of the two vectors
  mu.x <- mean(x)
  mu.y <- mean(y)

  # Calculate the variances of the two vectors
  var.x <- var(x)
  var.y <- var(y)

  # Calculate the standard deviations of the two vectors
  sd.x <- sqrt(var.x)
  sd.y <- sqrt(var.y)

  # Calculate the t statistic
  t.stat <- (mu.x - mu.y) / (sd.x / sqrt(length(x)) + sd.y / sqrt(length(y)))

  # Return the t statistic
  return(t.stat)
}

# Function to calculate the t statistic for a comparison of two means, assuming that the variances are equal
t.test.equal.var <- function(x, y) {

xbar <- mean(x)
ybar <- mean(y)
s1 <- sd(x)
s2 <- sd(y)

# Calculate the pooled standard deviation
s <- sqrt(((length(x) - 1) * s1^2 + (length(y) - 1) * s2^2) / (length(x) + length(y) - 2))

# Calculate the t statistic
t.stat <- (xbar - ybar) / (s * sqrt(1/length(x) + 1/length(y)))

  # Return the t statistic
  return(t.stat)
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

Results=data.frame(matrix(ncol=5,nrow=nrow(Input)))
colnames(Results)=c("DefTT","EqVarTT","MyTT","MyEqVarTT","LR")
UnderResults=Results
colnames(UnderResults)=sprintf("%sUnder",colnames(Results))
OverResults=Results
colnames(OverResults)=sprintf("%sOver",colnames(Results))
MeanResults=Results
colnames(MeanResults)=sprintf("%sMean",colnames(Results))

for (r in 1:nrow(Input)) {
	SimResults=data.frame(matrix(ncol=ncol(Results),nrow=Input$NSims[r]))
	colnames(SimResults)=colnames(Results)
	for (s in 1:Input$NSims[r]) {
		Seed=s+(r-1)*Input$NSims
		set.seed(Seed)
		Cases=rpois(Input$NCase[r],Input$CaseProb[r])
		Controls=rpois(Input$NCont[r],Input$ContProb[r])
		tt=t.test(Controls,Cases,var.equal=FALSE)
		p=tt$p.value
		SimResults$DefTT[s]=log10(p)*as.numeric(sign(tt$estimate[1]-tt$estimate[2]))
		tt=t.test(Controls,Cases,var.equal=TRUE)
		p=tt$p.value
		SimResults$EqVarTT[s]=log10(p)*as.numeric(sign(tt$estimate[1]-tt$estimate[2]))
		DF=Input$NCase[r]+Input$NCont[r]-2
		tt=t.test.unequal.var(Controls,Cases)
		p=2*pt(q=abs(tt), df=DF, lower.tail=FALSE) # one-tailed
		SimResults$MyTT[s]=log10(p)*sign(mean(Controls)-mean(Cases))
		tt=t.test.equal.var(Controls,Cases)
		p=2*pt(q=abs(tt), df=DF, lower.tail=FALSE)
		SimResults$MyEqVarTT[s]=log10(p)*sign(mean(Controls)-mean(Cases))
		ch2=LR.chisq(Controls,Cases)
		p=1-pchisq(q=ch2, df=1)
		SimResults$LR[s]=log10(p)*sign(mean(Controls)-mean(Cases))
	}
	write.table(SimResults,sprintf("%s.SimResults.%d.txt",TestName,r),row.names=FALSE,quote=FALSE,sep="\t")
	for (c in 1:ncol(SimResults)) {
		UnderResults[r,c]=length(which(SimResults[,c]<log10(0.05)))
		OverResults[r,c]=length(which(SimResults[,c]>-log10(0.05)))
		MeanResults[r,c]=mean(SimResults[,c])
	}
}

Results=cbind(UnderResults,OverResults)
Results=cbind(Results,MeanResults)


write.table(Results,ResultsFile,row.names=FALSE,quote=FALSE,sep="\t")


