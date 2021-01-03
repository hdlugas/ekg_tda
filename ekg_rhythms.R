###test for atrial fibrillation###
get_afib<-function(pr_int_lhe,rpeaks,npwaves){
  afib<-0
  count<-(npwaves-nrow(pr_int_lhe)) #count is the number of pr-intervals with no detected p-waves
  
  ##test 1
  if(npwaves<=(0.7*(nrow(rpeaks)-1))){
    afib<-1
  }
  
  ##test 2
  if(count>=(0.3*(nrow(rpeaks)-1))){
    afib<-1
  }
  
  ##test 3{
  if(npwaves<=(0.85*(nrow(rpeaks)-1)) && count>(0.15*(nrow(rpeaks)-1))){
    afib<-1
  }

  return(afib)
}



###test for atrial flutter###
get_atrial_flutter<-function(rpeaks,npwaves){
  atrial_flutter<-0
  if(npwaves>=(1.3*(nrow(rpeaks)-1))){
    atrial_flutter<-1
  }
  return(atrial_flutter)
}



###test for first-degree AV blockr###
get_first_avb<-function(npwaves,pr_interval_avg){
  first_avb<-0
  if(npwaves>(0.8*(nrow(rpeaks)-1)) && pr_interval_avg>0.2 && npwaves<(1.2*(nrow(rpeaks)-1))){
    first_avb<-1
  }
  return(first_avb)
}



###test for tachycardia###
get_tachycardia<-function(heart_rate){
  tachycardia<-0
  if(heart_rate>=100){
    tachycardia<-1
  }
  return(tachycardia)
}



###test for bradycardia###
get_bradycardia<-function(heart_rate){
  bradycardia<-0
  if(heart_rate<=60){
    bradycardia<-1
  }
  return(bradycardia)
}



#get center of mass of representative cycles given by m1 and m2
get_com<-function(m1, m2){
  com_coord<-c(0,0)
  xcoord1<-cbind(m1[,1])
  xcoord2<-cbind(m1[,2])
  xcoordt<-rbind(xcoord1,xcoord2)
  xcoord<-unique(xcoordt)
  ycoord1<-cbind(m2[,1])
  ycoord2<-cbind(m2[,2])
  ycoordt<-rbind(ycoord1,ycoord2)
  ycoord<-unique(ycoordt)
  com_coord[1]<-mean(xcoord)
  com_coord[2]<-mean(ycoord)
  return(com_coord)
}



#get the center of mass of each H1 feature with more than three data making up its 
#representative cycle along with the persistence of the feature
get_com_persistence<-function(diag){
  temp<-cbind(which(diag[[1]][,1]>0)) #this gives the indices of the rows of diag which correspond only to H1 and not H0 features
  persistence_H1_vec<-cbind()
  x_com_vec<-cbind()
  y_com_vec<-cbind()
  for(i in min(temp):n){
    if(nrow(diag[[4]][i][[1]][,,1])>=3){
      persistence_H1_vec<-append(persistence_H1_vec, diag[[1]][i,3]-diag[[1]][i,2], after=length(persistence_H1_vec))
      m1<-diag[[4]][i][[1]][,,1]
      m2<-diag[[4]][i][[1]][,,2]
      xcoord1<-cbind(m1[,1])
      xcoord2<-cbind(m1[,2])
      xcoordt<-rbind(xcoord1,xcoord2)
      xcoord<-unique(xcoordt)
      ycoord1<-cbind(m2[,1])
      ycoord2<-cbind(m2[,2])
      ycoordt<-rbind(ycoord1,ycoord2)
      ycoord<-unique(ycoordt)
      x_com_vec<-append(x_com_vec, mean(xcoord), after=length(x_com_vec))
      y_com_vec<-append(y_com_vec, mean(ycoord), after=length(y_com_vec))
    }
  }
  persistence_H1_vec<-cbind(persistence_H1_vec)
  
  com_persistence<-matrix(0,nrow(persistence_H1_vec),3)
  com_persistence[,1]<-x_com_vec
  com_persistence[,2]<-y_com_vec
  com_persistence[,3]<-persistence_H1_vec

  return(com_persistence)
}



#get all H1 features with a positive y center of mass coordinate within a given RR-interval
get_interval_H1_pos_y<-function(diag,com_persistence){
  interval_H1<-matrix(0,nrow(rpeaks)-1,20) 
  for(j in 2:nrow(rpeaks)){
    count<-1
    for(i in 1:nrow(com_persistence)){
      xc<-com_persistence[i,1]
      yc<-com_persistence[i,2]
      persistence<-com_persistence[i,3]
      if(yc>0 && xc>rpeaks[j-1,1] && xc<rpeaks[j,1]){
        interval_H1[j-1,count]<-persistence
        count<-count+1
      }
    }
  }
  return(interval_H1)
}



#measure the PR-intervals
get_pr_intervals<-function(pr_int_lhe,pr_int_rhe,rpeaks){
  pr_intervals<-cbind()
  for(j in 2:nrow(rpeaks)){
    for(i in 1:nrow(pr_int_lhe)){
      if(pr_int_lhe[i]>rpeaks[j-1,1] && pr_int_lhe[i]<rpeaks[j,1]){
        for(k in 1:nrow(pr_int_rhe)){
          if(pr_int_rhe[k]>rpeaks[j-1,1] && pr_int_rhe[k]<rpeaks[j,1]){
            diff<-pr_int_rhe[k]-pr_int_lhe[i]
            pr_intervals<-append(pr_intervals,diff,after=length(pr_intervals))
          }
        }
      }
    }
  }
  return(pr_intervals)
}




#measure the QRS-durations
get_qrs_dur<-function(qrs_dur_lhe,qrs_dur_rhe,rpeaks){
  qrs_dur<-matrix(0,(nrow(rpeaks)-2),1)
  for(j in 2:(nrow(rpeaks)-1)){
    qrs_dur[j-1]<-qrs_dur_rhe[j]-qrs_dur_lhe[j-1]
  }
  return(qrs_dur)
}



#measure the QT-intervals
get_qt_intervals<-function(qrs_dur_lhe,t_rhe_vec,rpeaks){
  qt_intervals<-cbind()
  for(j in 2:(nrow(rpeaks)-1)){
    for(i in 1:nrow(t_rhe_vec)){
      if(t_rhe_vec[i]>rpeaks[j,1] && t_rhe_vec[i]<rpeaks[j+1,1]){
        for(k in 1:nrow(qrs_dur_lhe)){
          if(qrs_dur_lhe[k]>rpeaks[j-1,1] && qrs_dur_lhe[k]<rpeaks[j,1]){
            diff<-t_rhe_vec[i]-qrs_dur_lhe[k]
            qt_intervals<-append(qt_intervals,diff,after=length(qt_intervals))
          }
        }
      }
    }
  }
  return(qt_intervals)
}



#measure the ST-intervals
get_st_intervals<-function(qrs_dur_rhe,t_rhe_vec,rpeaks){
  st_intervals<-cbind()
  for(j in 2:nrow(rpeaks)){
    for(i in 1:nrow(t_rhe_vec)){
      if(t_rhe_vec[i]>rpeaks[j-1,1] && t_rhe_vec[i]<rpeaks[j,1]){
        for(k in 1:nrow(qrs_dur_rhe)){
          if(qrs_dur_rhe[k]>rpeaks[j-1,1] && qrs_dur_rhe[k]<rpeaks[j,1]){
            diff<-t_rhe_vec[i]-qrs_dur_rhe[k]
            st_intervals<-append(st_intervals,diff,after=length(st_intervals))
          }
        }
      }
    }
  }
  return(st_intervals)
}

