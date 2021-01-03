###include baseline###
add_baseline<-function(ekgt){
  ekg<-matrix(0,nrow(ekgt),1)
  for(i in 1:nrow(ekgt)){
    if(i%%3==0){
      ekg[i]<-0
    }
    if(i%%3!=0){
      ekg[i]<-ekgt[i,1]
    }
  }
  return(ekg)
}



###normalize data###
normalize_data<-function(ekg){
  max_signal<-0.00
  for (i in 1:nrow(ekg)){
    if(ekg[i,1]>=max_signal){
      max_signal=ekg[i,1]
    }
  }
  
  nekg<-matrix(0,nrow(ekg),1)
  
  for (i in 1:nrow(ekg)){
    nekg[i,1]<-ekg[i,1]/max_signal
  }
  return(nekg)
}



###include time axis###
time_axis<-function(ekgt,nekg,sf){
  total_time<-nrow(ekgt)/sf
  t<-matrix(0,nrow(ekgt),1)
  
  for(i in 2:nrow(ekgt)){
    t[i,1]<-t[i-1,1]+total_time/nrow(ekgt)
  }
  fekg<-cbind(t,nekg)
  return(fekg)
}



###get heart rate###
get_heart_rate<-function(rpeaks){
  sum<-0
  for(i in 2:nrow(rpeaks)){
    sum<-sum+rpeaks[i,1]-rpeaks[i-1,1]
  }
  avg_rr_interval<-sum/(nrow(rpeaks)-1)
  heart_rate<-60/avg_rr_interval
  return(heart_rate)
}



###get maximum element N from vector x###
maxN <- function(x, N){
  len <- length(x)
  sort(x,partial=len-N+1)[len-N+1]
}



###define data to start and end with R-wave peaks###
data_trim<-function(ekg_temp,rpeaks){
  count1<-0
  count2<-0
  for(i in 1:nrow(ekg_temp)){
    if(ekg_temp[i,1]>=rpeaks[1,1] && ekg_temp[i,1]<=max(rpeaks[,1])){
      count1<-count1+1
    }
  }
  fekg<-matrix(0,count1,2)
  for(i in 1:nrow(ekg_temp)){
    if(ekg_temp[i,1]>rpeaks[1,1] && ekg_temp[i,1]<rpeaks[nrow(rpeaks),1]){
      count2<-count2+1
      fekg[count2,]<-ekg_temp[i,]
    }
  }
  for(i in 1:nrow(fekg)){
    if(fekg[i,1]==0 && fekg[i,2]==0){
      fekg[i,1]<-NA
      fekg[i,2]<-NA
    }
  }
  fekg<-na.omit(fekg)
  return(fekg)
}


