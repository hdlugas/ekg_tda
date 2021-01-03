###analyze T-waves###
analyze_t_waves<-function(p_com_vec,interval_H1_pos_y,diag,com_persistence,pthres_t_wave_lower,pthres_t_wave_upper,rpeaks){
  ntwaves<-0
  t_x_com_vec<-cbind()
  t_y_com_vec<-cbind()
  t_duration_vec<-cbind()
  t_rhe_vec<-cbind()
  for(j in 2:nrow(rpeaks)){
    for(i in 1:nrow(com_persistence)){
      xc<-com_persistence[i,1]
      yc<-com_persistence[i,2]
      persistence<-com_persistence[i,3]
      if(yc>0 && yc<0.4 && xc>rpeaks[j-1,1] && xc<(rpeaks[j,1]-0.25*(rpeaks[j,1]-rpeaks[j-1,1])) && 
         persistence==max(interval_H1_pos_y[j-1,]) | yc>0 && yc<0.4 && xc>rpeaks[j-1,1] && 
         xc<(rpeaks[j,1]-0.25*(rpeaks[j,1]-rpeaks[j-1,1])) &&  persistence>=pthres_t_wave_lower && 
         persistence<=pthres_t_wave_upper){
        ntwaves<-ntwaves+1
        points(xc,yc,type="p",pch=19,col="blue")
        t_x_com_vec<-append(t_x_com_vec,xc,after=length(t_x_com_vec))
        t_y_com_vec<-append(t_y_com_vec,yc,after=length(t_y_com_vec))
        for(k in 1:nrow(diag[[1]])){
          if(diag[[1]][k,1]==1){
            p<-diag[[1]][k,3]-diag[[1]][k,2]
            if(p==persistence){
              m1<-diag[[4]][k][[1]][,,1]
              m2<-diag[[4]][k][[1]][,,2]
              t_rhe_vec<-append(t_rhe_vec,max(m1),after=length(t_rhe_vec))
              t_duration<-max(m1)-min(m1)
              t_duration_vec<-append(t_duration_vec,t_duration,after=length(t_duration_vec))
              for (k in 1:nrow(m1)){
                segments(m1[k,1],m2[k,1],m1[k,2],m2[k,2],col="blue")
              }
            }
          }
        }
      }
    }
  }  
  t_com_vec<-cbind(t_x_com_vec,t_y_com_vec)
  avg_t_duration<-mean(t_duration_vec)
  sd_t_duration<-sd(t_duration_vec)
  t_info<-list(ntwaves,t_com_vec,avg_t_duration,sd_t_duration,t_rhe_vec)
  return(t_info)
}

