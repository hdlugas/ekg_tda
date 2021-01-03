library(TDA)
p_time<-cbind(1:100)
t_time<-cbind(1:76)
q_time<-cbind(1:20)
s_time<-cbind(1:30)
r_time<-cbind(1:20)
#time<-cbind(time/250)
pt<-cbind(time[1:100])
tt<-cbind(time[400:500])
p_width<-12  #width of gaussian p-wave at half-max
t_width<-10
p<-cbind(0.2*exp(-1*((p_time-50)^2)/(2*p_width*p_width)))
t<-cbind(0.4*exp(-1*((t_time-50)^2)/(2*t_width*t_width))) 
q1<-cbind(-0.003*q_time[1:10])
q2<-cbind(0.003*q_time[10:20]-nrow(q_time)/500-0.02)
q<-rbind(q1,q2)
s1<-cbind(-0.003*s_time[1:15])
s2<-cbind(0.003*s_time[15:30]-nrow(s_time)/500-0.03)
s<-rbind(s1,s2)
r1<-cbind(r_time[1:10]/10)
r2<-cbind(-1*r_time[1:10]/10+1)
r<-rbind(r1,r2)
ekgt<-rbind(p,q,r,s,t)
plot(ekgt)


#add baseline
ekg<-matrix(0,nrow(ekgt),1)
for(i in 1:nrow(ekgt)){
  if(i%%3==0){
    ekg[i]<-0
  }
  if(i%%3!=0){
    ekg[i]<-ekgt[i,1]
  }
}

time<-cbind(1:248)
n_time<-cbind(time/248)
ekg<-cbind(n_time,ekg)

plot(ekg,main="Processed EKG with r = 0.025",ylab="Normalized Voltage",xlab="Time (seconds)")
for(i in 1:nrow(ekg)){
  symbols(ekg[i,1],ekg[i,2],circles=0.025,inches=FALSE,add=TRUE)
}

plot(ekg,main="Geometric Cech Complex with r = 0.025",ylab="Normalized Voltage",xlab="Time (seconds)")
for(i in 1:nrow(ekg)){
  symbols(ekg[i,1],ekg[i,2],circles=0.025,inches=FALSE,add=TRUE,bg="green")
}

diag<-ripsDiag(ekg,1,0.025,dist="euclidean",library="Dionysus",location=TRUE) 
plot(diag[["diagram"]],main="Persistent Homology with r = 0.025")



x<-c(0,0,1,1)
y<-c(0,1,1,0)
z<-cbind(x,y)
plot(z,main="Geometric Cech Complex with r = 0.71")
for(i in 1:nrow(z)){
  symbols(z[i,1],z[i,2],circles=0.71,inches=FALSE,add=TRUE,bg="green")
}

diag<-ripsDiag(z,1,2,dist="euclidean",library="Dionysus",location=TRUE)
plot(diag[["diagram"]],main="Persistent Homology of Square")

