###analyze P-waves###
analyze_p_waves<-function(interval_H1_pos_y,diag,com_persistence,pthres_p_wave_lower,pthres_p_wave_upper,rpeaks){
  npwaves<-0
  p_x_com_vec<-cbind()
  p_y_com_vec<-cbind()
  p_duration_vec<-cbind()
  p_lhe_vec<-cbind()
  for(j in 2:nrow(rpeaks)){
    for(i in 1:nrow(com_persistence)){
      xc<-com_persistence[i,1]
      yc<-com_persistence[i,2]
      persistence<-com_persistence[i,3]
      if(yc>0 && xc>rpeaks[j-1,1] && xc<rpeaks[j,1] && persistence!=max(interval_H1_pos_y[j-1,])
         && persistence>=pthres_p_wave_lower && persistence<=pthres_p_wave_upper){
        npwaves<-npwaves+1
        points(xc,yc,type="p",pch=19,col="red")
        p_x_com_vec<-append(p_x_com_vec,xc,after=length(p_x_com_vec))
        p_y_com_vec<-append(p_y_com_vec,yc,after=length(p_y_com_vec))
        for(k in 1:nrow(diag[[1]])){
          if(diag[[1]][k,1]==1){
            p<-diag[[1]][k,3]-diag[[1]][k,2]
            if(p==persistence){
              m1<-diag[[4]][k][[1]][,,1]
              m2<-diag[[4]][k][[1]][,,2]
              p_lhe_vec<-append(p_lhe_vec,min(m1),after=length(p_lhe_vec))
              p_duration<-max(m1)-min(m1)
              p_duration_vec<-append(p_duration_vec,p_duration,after=length(p_duration_vec))
              for (k in 1:nrow(m1)){
                segments(m1[k,1],m2[k,1],m1[k,2],m2[k,2],col="red")
              }
            }
          }
        }
      }
    }
  }  
  p_com_vec<-cbind(p_x_com_vec,p_y_com_vec)
  avg_p_duration<-mean(p_duration_vec)
  sd_p_duration<-sd(p_duration_vec)
  p_info<-list(npwaves,p_com_vec,avg_p_duration,sd_p_duration,p_lhe_vec)
  return(p_info)
}


#get the leff-hand endpoints of the PR-intervals
get_pr_int_lhe<-function(p_lhe_vec,rpeaks){
  pr_int_lhe<-cbind()
  pr_int_lhe_total<-matrix(0,nrow(rpeaks)-1,5) #this will be defined such that the n-th row holds the lhe of all detected p-waves within the n-th RR-interval
  count<-1
  for(j in 2:nrow(rpeaks)){
    for(i in 1:nrow(p_lhe_vec)){
      if(p_lhe_vec[i]>rpeaks[j-1,1] && p_lhe_vec[i]<rpeaks[j,1]){
        pr_int_lhe_total[j-1,count]<-p_lhe_vec[i]
      }
    }
  }
  
  for(i in 1:nrow(pr_int_lhe_total)){
    if(max(pr_int_lhe_total[i,]) != 0)
      pr_int_lhe<-append(pr_int_lhe,max(pr_int_lhe_total[i,]),after=length(pr_int_lhe))
  }
 
  return(pr_int_lhe)
}


