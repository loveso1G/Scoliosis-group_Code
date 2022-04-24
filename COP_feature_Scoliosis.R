#
rm(list=ls())
gc()
# install.packages("pracma")#"#" must be removed for the first installation
# install.packages("EMD")#"#" must be removed for the first installation
# signal: COPx, COPy
t1=Sys.time();print(t1)
rawdata=read.csv("16975967_4c.csv",header=T);dim(rawdata);#Just change 16975967_4c and make changes according to the input data
tt=rawdata[,1];copx=rawdata[,2]; copy=rawdata[,3];
char2="_result.csv";
char1="16975967_4c";char3=paste(char1,char2, sep='');
# EMD
library(EMD)
E1=emd(copx,tt,boundary="wave");
E1.no=E1$nimf;
E1.imf=E1$imf
E1.residue=E1$residue

E2=emd(copy,tt,boundary="wave");
E2.no=E2$nimf;
E2.imf=E2$imf
E2.residue=E2$residue

# time domain features
Result_all=matrix(0,1,42);
	# Feature-1: RMSD
m0E1=sqrt(mean((copx-mean(copx))^2))            ;Result_all[1]=m0E1 #x.0.1
m1E1=sqrt(mean((E1.imf[,1]-mean(E1.imf[,1]))^2));Result_all[7]=m1E1 #x.1.1
m2E1=sqrt(mean((E1.imf[,2]-mean(E1.imf[,2]))^2));Result_all[13]=m2E1 #x.2.1
m3E1=sqrt(mean((E1.imf[,3]-mean(E1.imf[,3]))^2));Result_all[19]=m3E1 #x.3.1
m4E1=sqrt(mean((E1.imf[,4]-mean(E1.imf[,4]))^2));Result_all[25]=m4E1 #x.4.1
m5E1=sqrt(mean((E1.imf[,5]-mean(E1.imf[,5]))^2));Result_all[31]=m5E1 #x.5.1
m6E1=sqrt(mean((E1.imf[,6]-mean(E1.imf[,6]))^2));Result_all[37]=m6E1 #x.6.1

m0E2=sqrt(mean((copy-mean(copy))^2))            ;Result_all[2]=m0E2 #y.0.1
m1E2=sqrt(mean((E2.imf[,1]-mean(E2.imf[,1]))^2));Result_all[8]=m1E2#y.1.1
m2E2=sqrt(mean((E2.imf[,2]-mean(E2.imf[,2]))^2));Result_all[14]=m2E2#y.2.1
m3E2=sqrt(mean((E2.imf[,3]-mean(E2.imf[,3]))^2));Result_all[20]=m3E2#y.3.1
m4E2=sqrt(mean((E2.imf[,4]-mean(E2.imf[,4]))^2));Result_all[26]=m4E2#y.4.1
m5E2=sqrt(mean((E2.imf[,5]-mean(E2.imf[,5]))^2));Result_all[32]=m5E2#y.5.1
m6E2=sqrt(mean((E2.imf[,6]-mean(E2.imf[,6]))^2));Result_all[38]=m6E2#y.6.1

SD0E1=sd(copx);	  	
SD1E1=sd(E1.imf[,1]);
SD2E1=sd(E1.imf[,2]);
SD3E1=sd(E1.imf[,3]);
SD4E1=sd(E1.imf[,4]);
SD5E1=sd(E1.imf[,5]);
SD6E1=sd(E1.imf[,6]);

SD0E2=sd(copy)
SD1E2=sd(E2.imf[,1]);
SD2E2=sd(E2.imf[,2]);
SD3E2=sd(E2.imf[,3]);
SD4E2=sd(E2.imf[,4]);
SD5E2=sd(E2.imf[,5]);
SD6E2=sd(E2.imf[,6]);

# http://127.0.0.1:17734/library/pracma/html/entropy.html
	#Feature-2:  Approximate Entropy 
 library(pracma)
print("Approximate Entropy E1");Sys.time()
AP0E1=approx_entropy(copx,edim=2,r=0.2*SD0E1,elag = 1);	 Result_all[3]=AP0E1; #x.0.4
AP1E1=approx_entropy(E1.imf[,1],edim=2,r=0.2*SD1E1,elag = 1);Result_all[9]=AP1E1; #x.1.4
AP2E1=approx_entropy(E1.imf[,2],edim=2,r=0.2*SD2E1,elag = 1);Result_all[15]=AP2E1; #x.2.4
AP3E1=approx_entropy(E1.imf[,3],edim=2,r=0.2*SD3E1,elag = 1);Result_all[21]=AP3E1; #x.3.4
AP4E1=approx_entropy(E1.imf[,4],edim=2,r=0.2*SD4E1,elag = 1);Result_all[27]=AP4E1; #x.4.4
AP5E1=approx_entropy(E1.imf[,5],edim=2,r=0.2*SD5E1,elag = 1);Result_all[33]=AP5E1; #x.5.4
AP6E1=approx_entropy(E1.imf[,6],edim=2,r=0.2*SD6E1,elag = 1);Result_all[39]=AP6E1; #x.6.4


print("Approximate Entropy E2");Sys.time()
AP0E2=approx_entropy(copy,edim=2,r=0.2*SD0E2,elag = 1);	 Result_all[4]=AP0E2; #y.0.4
AP1E2=approx_entropy(E2.imf[,1],edim=2,r=0.2*SD1E2,elag = 1);Result_all[10]=AP1E2; #y.1.4
AP2E2=approx_entropy(E2.imf[,2],edim=2,r=0.2*SD2E2,elag = 1);Result_all[16]=AP2E2; #y.2.4
AP3E2=approx_entropy(E2.imf[,3],edim=2,r=0.2*SD3E2,elag = 1);Result_all[22]=AP3E2; #y.3.4
AP4E2=approx_entropy(E2.imf[,4],edim=2,r=0.2*SD4E2,elag = 1);Result_all[28]=AP4E2; #y.4.4
AP5E2=approx_entropy(E2.imf[,5],edim=2,r=0.2*SD5E2,elag = 1);Result_all[34]=AP5E2; #y.5.4
AP6E2=approx_entropy(E2.imf[,6],edim=2,r=0.2*SD6E2,elag = 1);Result_all[40]=AP6E2; #y.6.4


	#Feature-3:  Sample entropy

print("Sample Entropy E1");Sys.time()
SA0E1=sample_entropy(copx,edim=2,r=0.2*SD1E1,tau = 1);	Result_all[5]=SA0E1; #x.0.5
SA1E1=sample_entropy(E1.imf[,1],edim=2,r=0.2*SD1E1,tau = 1);Result_all[11]=SA1E1; #x.1.5
SA2E1=sample_entropy(E1.imf[,2],edim=2,r=0.2*SD2E1,tau = 1);Result_all[17]=SA2E1; #x.2.5
SA3E1=sample_entropy(E1.imf[,3],edim=2,r=0.2*SD3E1,tau = 1);Result_all[23]=SA3E1; #x.3.5
SA4E1=sample_entropy(E1.imf[,4],edim=2,r=0.2*SD4E1,tau = 1);Result_all[29]=SA4E1; #x.4.5
SA5E1=sample_entropy(E1.imf[,5],edim=2,r=0.2*SD5E1,tau = 1);Result_all[35]=SA5E1; #x.5.5
SA6E1=sample_entropy(E1.imf[,6],edim=2,r=0.2*SD6E1,tau = 1);Result_all[41]=SA6E1; #x.6.5

print("Sample Entropy E2");Sys.time()
SA0E2=sample_entropy(copy,edim=2,r=0.2*SD1E1,tau = 1);	Result_all[6]=SA0E2; #y.0.5
SA1E2=sample_entropy(E2.imf[,1],edim=2,r=0.2*SD1E2,tau = 1);Result_all[12]=SA1E2; #y.1.5
SA2E2=sample_entropy(E2.imf[,2],edim=2,r=0.2*SD2E2,tau = 1);Result_all[18]=SA2E2; #y.2.5
SA3E2=sample_entropy(E2.imf[,3],edim=2,r=0.2*SD3E2,tau = 1);Result_all[24]=SA3E2; #y.3.5
SA4E2=sample_entropy(E2.imf[,4],edim=2,r=0.2*SD4E2,tau = 1);Result_all[30]=SA4E2; #y.4.5
SA5E2=sample_entropy(E2.imf[,5],edim=2,r=0.2*SD5E2,tau = 1);Result_all[36]=SA5E2; #y.5.5
SA6E2=sample_entropy(E2.imf[,6],edim=2,r=0.2*SD6E2,tau = 1);Result_all[42]=SA6E2; #y.6.5




print(Result_all);
write.csv(Result_all,  file =char3)
print("Sample Entropy E2");Sys.time()