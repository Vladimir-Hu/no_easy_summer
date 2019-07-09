x0 <- (PA_CONST=10000,PR_CONST=10000,PA=1000,PR=1000,mRNA_A=0,mRNA_R=0,A=0,R=0)
spd <- c()
spd[1:88]=0
spd[45:46]=+1
spd[58:59]=+1
spd[71]=+1
spd[83]=+1
spd[52]=-1
spd[64]=-1
spd[73]=-1
spd[76]=-1
spd[88]=-1
a <- c("0.00022*{PA_CONST}","(0.1875*{A}**2)/(0.2**2+{A}**2)","0.000625*{PR_CONST}","(0.0625*{A}**5)/(0.12**5+{A}**5)","0.2*{mRNA_A}","0.2*{mRNA_R}","4e-8*{R}*{A}","0.005*{mRNA_A}","0.005*{mRNA_R}","0.0001*{A}","0.0001*{R}")
nu <- matrix(spd,nrow=8,byrow=T)
outssa <- ssa(x0,a,nu,tf=100,method="D")
outtau <- ssa(x0,a,nu,tf=100,method="ETL")
outotl <- ssa(x0,a,nu,tf=100,method="OTL")