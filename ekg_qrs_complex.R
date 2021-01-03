###get number of r peaks###
get_number_r_peaks<-function(ekg_temp){
  nrpeaks<-0
  for(i in 4:nrow(ekg_temp)-3){
    if(ekg_temp[i,2]>=0.5 && ekg_temp[i,2]>=ekg_temp[i+3,2] && ekg_temp[i,2]>=ekg_temp[i-3,2] && 
       ekg_temp[i,2]>=ekg_temp[i-2,2] && ekg_temp[i,2]>=ekg_temp[i+2,2] && ekg_temp[i,2]>=ekg_temp[i+1,2] && ekg_temp[i,2]>=ekg_temp[i-1,2]){
      nrpeaks<-nrpeaks+1
    }
  }
  return(nrpeaks)
}



###get location of R-wave peaks###
get_rpeak_coord<-function(ekg_temp,nrpeaks){
  count<-0
  rpeaks<-matrix(0,nrpeaks,2)
  for(i in 4:nrow(ekg_temp)-3){
    if(ekg_temp[i,2]>=0.5 && ekg_temp[i,2]>=ekg_temp[i+3,2] && ekg_temp[i,2]>=ekg_temp[i-3,2] && 
       ekg_temp[i,2]>=ekg_temp[i-2,2] && ekg_temp[i,2]>=ekg_temp[i+2,2] && ekg_temp[i,2]>=ekg_temp[i+1,2] && ekg_temp[i,2]>=ekg_temp[i-1,2]){
      count<-count+1
      rpeaks[count,]<-ekg_temp[i,]
    }
  }
  return(rpeaks)
}



###get RR-intervals###
get_rr_intervals<-function(rpeaks){
  rr_intervals<-matrix(0,nrow(rpeaks)-1,1)
  for(i in 2:nrow(rpeaks)){
    rr_intervals[i-1]<-rpeaks[i,1]-rpeaks[i-1,1]
  }
  return(rr_intervals)
}



###get the standard deviation of the RR-intervals and then propagate the error towards the heart rate###
get_heart_rate_error<-function(heart_rate,rr_intervals){
  heart_rate_error<-heart_rate*sd(rr_intervals)/mean(rr_intervals)
  if(heart_rate_error<=0.50){
    heart_rate_error<-1
  }
  heart_rate_error<-sprintf("%.0f",heart_rate_error)
  return(heart_rate_error)
}



#analyze q-waves
analyze_q_waves<-function(interval_H1_pos_y,diag,com_persistence,pthres_q_wave_lower,pthres_q_wave_upper,rpeaks){
  #first get all H1 features with persistence and COM within ranges for q-waves
  q_x_com_vec_raw<-matrix(0,(nrow(rpeaks)-1),5)
  q_y_com_vec_raw<-matrix(0,(nrow(rpeaks)-1),5)
  q_lhe_vec_raw<-matrix(0,(nrow(rpeaks)-1),5)
  for(j in 2:nrow(rpeaks)){
    count<-1
    for(i in 1:nrow(com_persistence)){
      xc<-com_persistence[i,1]
      yc<-com_persistence[i,2]
      persistence<-com_persistence[i,3]
      if(yc<0 && xc>(rpeaks[j-1,1]+0.75*(rpeaks[j,1]-rpeaks[j-1,1])) && xc<rpeaks[j,1]
         && persistence>=pthres_q_wave_lower && persistence<=pthres_q_wave_upper){
        q_x_com_vec_raw[j-1,count]<-xc
        q_y_com_vec_raw[j-1,count]<-yc
        count<-count+1
      }
    }
  }  
  
  #now filter the features detected above by considering only the right-most  
  #feature meeting the criteria within an RR-interval to be a Q-wave
  nqwaves<-0
  q_x_com_vec<-cbind()
  q_y_com_vec<-cbind()
  q_lhe_vec<-cbind()
  for(j in 2:nrow(rpeaks)){
    for(i in 1:nrow(com_persistence)){
      xc<-com_persistence[i,1]
      yc<-com_persistence[i,2]
      persistence<-com_persistence[i,3]
      if(xc==max(q_x_com_vec_raw[j-1])){
        nqwaves<-nqwaves+1
        points(xc,yc,type="p",pch=19,col="green")
        q_x_com_vec<-append(q_x_com_vec,xc,after=length(q_x_com_vec))
        q_y_com_vec<-append(q_y_com_vec,yc,after=length(q_y_com_vec))
        for(k in 1:nrow(diag[[1]])){
          if(diag[[1]][k,1]==1){
            p<-diag[[1]][k,3]-diag[[1]][k,2]
            if(p==persistence){
              m1<-diag[[4]][k][[1]][,,1]
              m2<-diag[[4]][k][[1]][,,2]
              q_lhe_vec<-append(q_lhe_vec,min(m1),after=length(q_lhe_vec))
              for (k in 1:nrow(m1)){
                segments(m1[k,1],m2[k,1],m1[k,2],m2[k,2],col="green")
              }
            }
          }
        }
      }
    }
  }  
  q_com_vec<-cbind(q_x_com_vec,q_y_com_vec)
  q_info<-list(nqwaves,q_com_vec,q_lhe_vec)
  return(q_info)
}




#get the rhe's of the PR-interval/lhe's of the QT-interval
get_pr_int_rhe<-function(q_lhe_vec,rpeaks,ekg_temp){
  pr_int_rhe<-matrix(0,(nrow(rpeaks)-1),1)
  
  #first populate pr_int_rhe with x-coordinate of the right-most point to the left of a QRS-complex under isoelectric line
  for(j in 2:nrow(rpeaks)){
    index<-which(ekg_temp[,1]==rpeaks[j,1])
    if(ekg_temp[index-1,2]<0){
      pr_int_rhe[j-1]<-ekg_temp[index-1,1]
    }
    else if(ekg_temp[index-2,2]<0){
      pr_int_rhe[j-1]<-ekg_temp[index-2,1]
    }
    else if(ekg_temp[index-3,2]<0){
      pr_int_rhe[j-1]<-ekg_temp[index-3,1]
    }
    else if(ekg_temp[index-4,2]<0){
      pr_int_rhe[j-1]<-ekg_temp[index-4,1]
    }
    else if(ekg_temp[index-5,2]<0){
      pr_int_rhe[j-1]<-ekg_temp[index-5,1]
    }
    else if(ekg_temp[index-6,2]<0){
      pr_int_rhe[j-1]<-ekg_temp[index-6,1]
    }
    else if(ekg_temp[index-7,2]<0){
      pr_int_rhe[j-1]<-ekg_temp[index-7,1]
    }
    else if(ekg_temp[index-8,2]<0){
      pr_int_rhe[j-1]<-ekg_temp[index-8,1]
    }
    else if(ekg_temp[index-9,2]<0){
      pr_int_rhe[j-1]<-ekg_temp[index-9,1]
    }
    else if(ekg_temp[index-10,2]<0){
      pr_int_rhe[j-1]<-ekg_temp[index-10,1]
    }
    else if(ekg_temp[index-11,2]<0){
      pr_int_rhe[j-1]<-ekg_temp[index-11,1]
    }
    else if(ekg_temp[index-12,2]<0){
      pr_int_rhe[j-1]<-ekg_temp[index-12,1]
    }
  }
  
  #now overwrite pr_int_rhe for detected Q-waves 
  for(j in 2:nrow(rpeaks)){
    for(i in 1:nrow(q_lhe_vec)){
      if(q_lhe_vec[i]>rpeaks[j-1,1] && q_lhe_vec[i]<rpeaks[j,1]){
        pr_int_rhe[j-1]<-q_lhe_vec[i]
      }
    }
  }

  return(pr_int_rhe)
}




#analyze s-waves
analyze_s_waves<-function(interval_H1_pos_y,diag,com_persistence,pthres_s_wave_lower,pthres_s_wave_upper,rpeaks){
  #first get all H1 features with persistence and COM within ranges for s-waves
  s_x_com_vec_raw<-matrix(0,(nrow(rpeaks)-1),5)
  s_y_com_vec_raw<-matrix(0,(nrow(rpeaks)-1),5)
  for(j in 2:nrow(rpeaks)){
    count<-1
    for(i in 1:nrow(com_persistence)){
      xc<-com_persistence[i,1]
      yc<-com_persistence[i,2]
      persistence<-com_persistence[i,3]
      if(yc<0 && xc>rpeaks[j-1,1] && xc<(rpeaks[j,1]-0.75*(rpeaks[j,1]-rpeaks[j-1,1]))
         && persistence>=pthres_s_wave_lower && persistence<=pthres_s_wave_upper){
        s_x_com_vec_raw[j-1,count]<-xc
        s_y_com_vec_raw[j-1,count]<-yc
        count<-count+1
      }
    }
  }  
  
  #now filter the features detected above by considering only the right-most  
  #feature meeting the criteria within an RR-interval to be a Q-wave
  nswaves<-0
  s_x_com_vec<-cbind()
  s_y_com_vec<-cbind()
  s_rhe_vec<-cbind()
  for(j in 2:nrow(rpeaks)){
    for(i in 1:nrow(com_persistence)){
      xc<-com_persistence[i,1]
      yc<-com_persistence[i,2]
      persistence<-com_persistence[i,3]
      if(xc==max(s_x_com_vec_raw[j-1]) && yc==max(s_y_com_vec_raw[j-1])){
        nswaves<-nswaves+1
        points(xc,yc,type="p",pch=19,col="purple")
        s_x_com_vec<-append(s_x_com_vec,xc,after=length(s_x_com_vec))
        s_y_com_vec<-append(s_y_com_vec,yc,after=length(s_y_com_vec))
        for(k in 1:nrow(diag[[1]])){
          if(diag[[1]][k,1]==1){
            p<-diag[[1]][k,3]-diag[[1]][k,2]
            if(p==persistence){
              m1<-diag[[4]][k][[1]][,,1]
              m2<-diag[[4]][k][[1]][,,2]
              s_rhe_vec<-append(s_rhe_vec,max(m1),after=length(s_rhe_vec))
              for (k in 1:nrow(m1)){
                segments(m1[k,1],m2[k,1],m1[k,2],m2[k,2],col="purple")
              }
            }
          }
        }
      }
    }
  }  
  s_com_vec<-cbind(s_x_com_vec,s_y_com_vec)
  s_info<-list(nswaves,s_com_vec,s_rhe_vec)
  return(s_info)
}




#get the rhe's of the QRS-duration/lhe's of the ST-interval
get_st_int_lhe<-function(s_rhe_vec,rpeaks,ekg_temp){
  st_int_lhe<-matrix(0,(nrow(rpeaks)-1),1)
  
  #first populate st_int_lhe with x-coordinate of the left-most point to the right of a QRS-complex under isoelectric line
  for(j in 2:nrow(rpeaks)){
    index<-which(ekg_temp[,1]==rpeaks[j,1])
    if(ekg_temp[index+1,2]<0){
      st_int_lhe[j-1]<-ekg_temp[index-1,1]
    }
    else if(ekg_temp[index+2,2]<0){
      st_int_lhe[j-1]<-ekg_temp[index-2,1]
    }
    else if(ekg_temp[index+3,2]<0){
      st_int_lhe[j-1]<-ekg_temp[index-3,1]
    }
    else if(ekg_temp[index+4,2]<0){
      st_int_lhe[j-1]<-ekg_temp[index-4,1]
    }
    else if(ekg_temp[index+5,2]<0){
      st_int_lhe[j-1]<-ekg_temp[index-5,1]
    }
    else if(ekg_temp[index+6,2]<0){
      st_int_lhe[j-1]<-ekg_temp[index-6,1]
    }
    else if(ekg_temp[index+7,2]<0){
      st_int_lhe[j-1]<-ekg_temp[index-7,1]
    }
    else if(ekg_temp[index+8,2]<0){
      st_int_lhe[j-1]<-ekg_temp[index-8,1]
    }
    else if(ekg_temp[index+9,2]<0){
      st_int_lhe[j-1]<-ekg_temp[index-9,1]
    }
    else if(ekg_temp[index+10,2]<0){
      st_int_lhe[j-1]<-ekg_temp[index-10,1]
    }
    else if(ekg_temp[index+11,2]<0){
      st_int_lhe[j-1]<-ekg_temp[index-11,1]
    }
    else if(ekg_temp[index+12,2]<0){
      st_int_lhe[j-1]<-ekg_temp[index-12,1]
    }
  }
  
  #now overwrite st_int_lhe for detected S-waves 
  for(j in 2:nrow(rpeaks)){
    for(i in 1:nrow(s_rhe_vec)){
      if(s_rhe_vec[i]>rpeaks[j-1,1] && s_rhe_vec[i]<rpeaks[j,1]){
        st_int_lhe[j-1]<-s_rhe_vec[i]
      }
    }
  }
  
  return(st_int_lhe)
}


