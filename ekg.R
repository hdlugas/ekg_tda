#author: Hunter Dlugas
#email: hunter.dlugas@ucsf.edu
#This is an algorithm which detects EKG features using persistent homology and then measures
#the PR,QT,ST-interval, QRS-duration, P-wave duration, and T-wave duration based off of the
#upper and lower time axis bounds of the representative cycles of the detected features within 
#the dataset. There are also flags which indicate atrial fibrillation, atrial flutter, and 
#first-degree AV block


library(TDA)
source("C:/Users/Hunter/Downloads/EKG_Project/ekg_preprocessing_functions.R")
source("C:/Users/Hunter/Downloads/EKG_Project/ekg_p_waves.R")
source("C:/Users/Hunter/Downloads/EKG_Project/ekg_qrs_complex.R")
source("C:/Users/Hunter/Downloads/EKG_Project/ekg_t_waves.R")
source("C:/Users/Hunter/Downloads/EKG_Project/ekg_rhythms.R")
setwd("C:/Users/Hunter/Downloads/EKG_Project")
#raw_ekg<-read.csv('afib_ekg_data1.csv')
# raw_ekg<-read.csv('afib_ekg_data2.csv')
# raw_ekg<-read.csv('afib_ekg_data3.csv')
# raw_ekg<-read.csv('afib_ekg_data4.csv')
# raw_ekg<-read.csv('afib_ekg_data5.csv')
# raw_ekg<-read.csv('afib_ekg_data6.csv')
# raw_ekg<-read.csv('afib_ekg_data7.csv')
# raw_ekg<-read.csv('afib_ekg_data8.csv')
# raw_ekg<-read.csv('afib_ekg_data9.csv')
# raw_ekg<-read.csv('afib_ekg_data10.csv')
# raw_ekg<-read.csv('af_ekg_data1.csv')
raw_ekg<-read.csv('sr_ekg_data1.csv')


###user input###
pthres_p_wave_lower<-0.017
pthres_p_wave_upper<-0.07
pthres_q_wave_lower<-0.001
pthres_q_wave_upper<-0.05
pthres_s_wave_lower<-0.001
pthres_s_wave_upper<-0.06
pthres_t_wave_lower<-0.07
pthres_t_wave_upper<-0.1
sf<-500
#note that a neural network may be constructed to better determine these thress. this may be difficult
#since it would likely require identification of each wave by eye to generate training data


###preprocessing of data filtered from python script###
ekgt<-cbind(raw_ekg[1:nrow(raw_ekg),2])
ekg<-add_baseline(ekgt)
nekg<-normalize_data(ekg)
ekg_temp<-time_axis(ekgt,nekg,sf)
nrpeaks<-get_number_r_peaks(ekg_temp)
rpeaks<-get_rpeak_coord(ekg_temp,nrpeaks)
fekg<-data_trim(ekg_temp,rpeaks)
heart_rate<-get_heart_rate(rpeaks)
rr_intervals<-get_rr_intervals(rpeaks)
heart_rate_error<-get_heart_rate_error(heart_rate,rr_intervals)



###compute homology###
diag<-ripsDiag(fekg,1,0.14,dist="euclidean",library="Dionysus",location=TRUE)
com_persistence<-get_com_persistence(diag) #this is an nx3 matrix with the n being the number of detected 
#H1 features.  The first two columns are the x and y coordinate of the center of mass of the representative
#cycle of the H1 feature, and the third column is the persistence of that H1 feature
interval_H1_pos_y<-get_interval_H1_pos_y(diag,com_persistence) ##note that the n-th row of interval_H1 will be defined to 
#hold the persistence values of the H1 features with yc>0 found in the n-th R-R interval




###make plots###
#png(filename="/wsu/home/fy/fy73/fy7392/ekg_homology.png")
plot(diag[["diagram"]],main="Patient One: Persistent Homology")
#dev.off()

#png(filename="/wsu/home/fy/fy73/fy7392/ekg_data.png")
plot(fekg,main="Patient One: Processed Data with Cycle Representations",xlab="Time (seconds)",ylab="Normalized Voltage (arbitrary units)")
points(rpeaks,pch=19,col="green")
#dev.off()




###analyze p-waves###
p_info<-analyze_p_waves(interval_H1_pos_y,diag,com_persistence,pthres_p_wave_lower,pthres_p_wave_upper,rpeaks)
npwaves<-p_info[1][[1]]
p_com_vec<-cbind(p_info[2][[1]])
avg_p_duration<-p_info[3][[1]]
sd_p_duration<-p_info[4][[1]]
p_lhe_vec<-cbind(p_info[5][[1]])
pr_int_lhe<-cbind(get_pr_int_lhe(p_lhe_vec,rpeaks)) #this is a vector of the lhe's of pr-intervals. There is
#at most one detected lhe within each rr-interval.  There may also be none.




###analyze q-waves###
q_info<-analyze_q_waves(interval_H1_pos_y,diag,com_persistence,pthres_q_wave_lower,pthres_q_wave_upper,rpeaks)
nqwaves<-q_info[1][[1]]
q_com_vec<-cbind(q_info[2][[1]])
q_lhe_vec<-cbind(q_info[3][[1]])
pr_int_rhe<-cbind(get_pr_int_rhe(q_lhe_vec,rpeaks,ekg_temp))




###analyze s-waves###
s_info<-analyze_s_waves(interval_H1_pos_y,diag,com_persistence,pthres_s_wave_lower,pthres_s_wave_upper,rpeaks)
nswaves<-s_info[1][[1]]
s_com_vec<-cbind(s_info[2][[1]])
s_rhe_vec<-cbind(s_info[3][[1]])
st_int_rhe<-cbind(get_st_int_lhe(s_rhe_vec,rpeaks,ekg_temp))




###analyze t-waves###
t_info<-analyze_t_waves(p_com_vec,interval_H1_pos_y,diag,com_persistence,pthres_t_wave_lower,pthres_t_wave_upper,rpeaks)
ntwaves<-t_info[1][[1]]
t_com_vec<-cbind(t_info[2][[1]])
avg_t_duration<-t_info[3][[1]]
sd_t_duration<-t_info[4][[1]]
t_rhe_vec<-cbind(t_info[5][[1]])




###calculate intervals###
pr_intervals<-cbind(get_pr_intervals(pr_int_lhe,pr_int_rhe,rpeaks))
pr_interval_avg<-mean(pr_intervals)
pr_interval_sd<-sd(pr_intervals)

qrs_dur_lhe<-pr_int_rhe
qrs_dur_rhe<-st_int_rhe
qrs_durations<-cbind(get_qrs_dur(qrs_dur_lhe,qrs_dur_rhe,rpeaks))
qrs_dur_avg<-mean(qrs_durations)
qrs_dur_sd<-sd(qrs_durations)

qt_intervals<-cbind(get_qt_intervals(qrs_dur_lhe,t_rhe_vec,rpeaks))
qt_interval_avg<-mean(qt_intervals)
qt_interval_sd<-sd(qt_intervals)

st_intervals<-cbind(get_st_intervals(qrs_dur_rhe,t_rhe_vec,rpeaks))
st_interval_avg<-mean(st_intervals)
st_interval_sd<-sd(st_intervals)




###classify rhythm###
#note that rhythm=0 if rhythm is not present and rhythm=1 if rhythm is present
afib<-get_afib(pr_int_lhe,rpeaks,npwaves)
atrial_flutter<-get_atrial_flutter(rpeaks,npwaves)
first_avb<-get_first_avb(npwaves,pr_interval_avg)
tachycardia<-get_tachycardia(heart_rate)
bradycardia<-get_bradycardia(heart_rate)




###write output to text file###
outfile<-"Heart_Rhythm_Summary.txt"
cat("Heart Rate ± SD: (",sprintf("%.0f",heart_rate),"±",heart_rate_error,")bpm","\n",file=outfile,append=TRUE)
cat("Number of P-waves detected:",npwaves,"\n",file=outfile,append=TRUE)
cat("Number of Q-waves detected:",nqwaves,"\n",file=outfile,append=TRUE)
cat("Number of R-waves detected:",nrow(rpeaks),"\n",file=outfile,append=TRUE)
cat("Number of S-waves detected:",nswaves,"\n",file=outfile,append=TRUE)
cat("Number of T-waves detected:",ntwaves,"\n",file=outfile,append=TRUE)
cat("Average PR-Interval ± SD: (",sprintf("%.3f",pr_interval_avg),"±",sprintf("%.3f",pr_interval_sd),")s","\n",file=outfile,append=TRUE)
cat("Average QRS-Duration ± SD: (",sprintf("%.3f",qrs_dur_avg),"±",sprintf("%.3f",qrs_dur_sd),")s","\n",file=outfile,append=TRUE)
cat("Average QT-Interval ± SD: (",sprintf("%.3f",qt_interval_avg),"±",sprintf("%.3f",qt_interval_sd),")s","\n",file=outfile,append=TRUE)
cat("Average ST-Interval ± SD: (",sprintf("%.3f",st_interval_avg),"±",sprintf("%.3f",st_interval_sd),")s","\n",file=outfile,append=TRUE)
cat("Average P-wave duration ± SD: (",sprintf("%.3f",avg_p_duration),"±",sprintf("%.3f",sd_p_duration),")s","\n",file=outfile,append=TRUE)
cat("Average T-wave duration ± SD: (",sprintf("%.3f",avg_t_duration),"±",sprintf("%.3f",sd_t_duration),")s","\n",file=outfile,append=TRUE)
cat("Atrial Fibrillation (0 is no and 1 is yes):",afib,"\n",file=outfile,append=TRUE)
cat("Atrial Flutter:",atrial_flutter,"\n",file=outfile,append=TRUE)
cat("First-Degree Atrioventricular Block:",first_avb,"\n",file=outfile,append=TRUE)
cat("Bradycardia:",bradycardia,"\n",file=outfile,append=TRUE)
cat("Tachycardia:",tachycardia,"\n",file=outfile,append=TRUE)

