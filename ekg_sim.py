#This script identifies P,Q,S, and T-waves of simulated ECG signals as optimal 1-cycles with certain properties, if 1-cycles with these properties exist.
#The PR-interval, QT-interval, ST-segment, QRS-duration, P-wave duration, and T-wave duration are then measured using these 1-cycles identified as P,Q,S, and T-waves
#and compared to the interval measurements determined by the ECG simulator.

import numpy as np
import homcloud.interface as hc
from numpy.random.mtrand import normal
import pandas as pd
from scipy import signal
from matplotlib.backends.backend_pdf import PdfPages
from processing import *
from cycles import *
from intervals import *

#input simulated ECG signal and interval measurements determined by ECG simulator
data=np.asarray(pd.read_csv('/home/hunter/ekg/data/ekg_sim_data.csv',delimiter=','))
intervals=np.asarray(pd.read_csv('/home/hunter/ekg/data/ekg_sim_intervals.csv',delimiter=','))#order: pr,qt,st,qrs,p,t

#create noise vector
#noise=np.random.normal(0,0.008,len(data[0,:]))

sf=500                                 #sampling frequency
n_samples=2                            #number of simulated ECG signals to analyze; max of 2000 for data provided on GitHub repo
x=np.linspace(0,4.0,2*len(data[0,:]))  #array of time coordinates of ECG signal


#declare empty arrays of interval measurements determined by ECG simulator
pr_sim=[]
qt_sim=[]
st_sim=[]
qrs_sim=[]
p_sim=[]
t_sim=[]
int_tda=np.zeros((len(data[:,0]),12))


#create pdf for images of optimal cycles drawn on ECG signals
pp = PdfPages('ekg_simulation_plots.pdf')
for i in range(0,n_samples):
    print(i)

    #import the i-th simulated ECG signal
    data_temp=data[i,:]

    #double the point density of the ECG signal
    data_temp=add_points(data_temp)

    #make ECG signal a nx2 matrix rather than nx1 array with specified sampling frequency to allow for the computation of 1-dimensional homological features
    ekg=np.zeros((len(data_temp),2))
    ekg[:,0]=x[:]
    ekg[:,1]=data_temp[:]

    #normalize the signal
    ekg[:,1]=normalize(ekg[:,1]) 

    #include an isoelectric baseline
    ekg=add_isoelectric_line(ekg)

    #compute R-wave peak coordinates
    rpeaks=signal.find_peaks(ekg[:,1],height=0.5,distance=10)
    r_peak_xcs,r_peak_xc_idx=get_rpeak_xcs(rpeaks,ekg)
    rr_int_avg,rr_int_sd=get_rr_interval(r_peak_xcs)

    #trim signal to begin and end with an R-wave peak
    ekg=trim(ekg,r_peak_xcs)

    #compute persistent homology of processed signal
    output=hc.PDList.from_alpha_filtration(ekg,no_squared=True,save_boundary_map=True,save_phtrees=True,save_to="pointcloud.pdgm")
    pd1=hc.PDList("pointcloud.pdgm").dth_diagram(1)
    persist=np.asarray(pd1.deaths-pd1.births)
    births=np.asarray(pd1.births)
    deaths=np.asarray(pd1.deaths)

    #compute optimal cycle representatives
    #cycle_xcs,cycle_ycs,bps=get_vol_opt_cycle_centroid_coords(persist,pd1)
    cycle_xcs,cycle_ycs,bps=get_card_opt_cycle_centroid_coords(persist,pd1)
    #cycle_xcs,cycle_ycs,bps=get_stab_vol_cycle_centroid_coords(persist,pd1)
    #cycle_xcs,cycle_ycs,bps=get_stab_subvol_cycle_centroid_coords(persist,pd1)

    #measure intervals of interest and compute the index locations of P,Q,S, and T-waves within the 1d array persist
    int_tda[i,:],idx_p,idx_q,idx_s,idx_t=get_intervals_and_H1wave_idxs(ekg,r_peak_xcs,rr_int_avg,persist,births,cycle_xcs,cycle_ycs,bps)
    #get_intervals_draw_cycles(ekg,r_peak_xcs,rr_int_avg,persist,births,stab_vol_cycle_xcs,stab_vol_cycle_ycs,bps)
    #get_intervals_draw_cycles(ekg,r_peak_xcs,rr_int_avg,persist,births,stab_subvol_cycle_xcs,stab_subvol_cycle_ycs,bps)

    #append the measurement of each interval to a vector
    if np.isnan(int_tda[i,0])==False:
        pr_sim.append(intervals[i,0])
    if np.isnan(int_tda[i,2])==False:
        qt_sim.append(intervals[i,1])
    if np.isnan(int_tda[i,4])==False:
        st_sim.append(intervals[i,2])
    if np.isnan(int_tda[i,6])==False:
        qrs_sim.append(intervals[i,3])
    if np.isnan(int_tda[i,8])==False:
        p_sim.append(intervals[i,4])
    if np.isnan(int_tda[i,10])==False:
        t_sim.append(intervals[i,5])
    
    #draw representative cycles
    count=i+1
    fig=draw_cycles(ekg,idx_p,idx_q,idx_s,idx_t,r_peak_xcs,cycle_xcs,bps,count)
    pp.savefig(fig)
pp.close()


#declare empty arrays of interval measurements using TDA
pr_tda=[]
qt_tda=[]
st_tda=[]
qrs_tda=[]
p_tda=[]
t_tda=[]

for i in range(0,n_samples):
    if np.isnan(int_tda[i,0])==False:
        pr_tda.append(int_tda[i,0])
    if np.isnan(int_tda[i,2])==False:
        qt_tda.append(int_tda[i,2])
    if np.isnan(int_tda[i,4])==False:
        st_tda.append(int_tda[i,4])
    if np.isnan(int_tda[i,6])==False:
        qrs_tda.append(int_tda[i,6])
    if np.isnan(int_tda[i,8])==False:
        p_tda.append(int_tda[i,8])
    if np.isnan(int_tda[i,10])==False:
        t_tda.append(int_tda[i,10])


#compute percent difference between interval measurements performed using TDA and the ECG simulator's parameters
pr_error_vec=[]
qt_error_vec=[]
st_error_vec=[]
qrs_error_vec=[]
p_error_vec=[]
t_error_vec=[]

for i in range(0,len(pr_tda)):
    pr_error_vec.append(np.abs(pr_tda[i]-pr_sim[i])/pr_sim[i]*100)

for i in range(0,len(qt_tda)):
    qt_error_vec.append(np.abs(qt_tda[i]-qt_sim[i])/qt_sim[i]*100)

for i in range(0,len(st_tda)):
    st_error_vec.append(np.abs(st_tda[i]-st_sim[i])/st_sim[i]*100)

for i in range(0,len(qrs_tda)):
    qrs_error_vec.append(np.abs(qrs_tda[i]-qrs_sim[i])/qrs_sim[i]*100)

for i in range(0,len(p_tda)):
    p_error_vec.append(np.abs(p_tda[i]-p_sim[i])/p_sim[i]*100)

for i in range(0,len(t_tda)):
    t_error_vec.append(np.abs(t_tda[i]-t_sim[i])/t_sim[i]*100)


#create a histogram for each interval of interest which shows the distribution of percent difference values for each measurement of the given 
#interval of interest performed using TDA and the ECG simulator's parameters
pp = PdfPages('error_histogram.pdf')
fig = plt.figure(figsize=[15,14])
ax1 = fig.add_subplot(2, 3, 1)
ax2 = fig.add_subplot(2, 3, 2)
ax3 = fig.add_subplot(2, 3, 3)
ax4 = fig.add_subplot(2, 3, 4)
ax5 = fig.add_subplot(2, 3, 5)
ax6 = fig.add_subplot(2, 3, 6)

ax1.hist(pr_error_vec,bins=np.linspace(0,100,21))
ax1.set_xlabel('Percent Difference (%)')
ax1.set_ylabel('Count')
ax1.set_title('Percent Difference Distribution\n of PR-Interval Measurements')

ax2.hist(qt_error_vec,bins=np.linspace(0,100,21))
ax2.set_xlabel('Percent Difference (%)')
ax2.set_ylabel('Count')
ax2.set_title('Percent Difference Distribution\n of QT-Interval Measurements')

ax3.hist(st_error_vec,bins=np.linspace(0,100,21))
ax3.set_xlabel('Percent Difference (%)')
ax3.set_ylabel('Count')
ax3.set_title('Percent Difference Distribution\n of ST-Segment Measurements')

ax4.hist(qrs_error_vec,bins=np.linspace(0,100,21))
ax4.set_xlabel('Percent Difference (%)')
ax4.set_ylabel('Count')
ax4.set_title('Percent Difference Distribution\n of QRS-Duration Measurements')

ax5.hist(p_error_vec,bins=np.linspace(0,100,21))
ax5.set_xlabel('Percent Difference (%)')
ax5.set_ylabel('Count')
ax5.set_title('Percent Difference Distribution\n of P-wave Duration Measurements')

ax6.hist(t_error_vec,bins=np.linspace(0,100,21))
ax6.set_xlabel('Percent Difference (%)')
ax6.set_ylabel('Count')
ax6.set_title('Percent Difference Distribution\n of T-wave Duration Measurements')

pp.savefig(fig)
pp.close()


#create a text file with the mean and standard deviation of the percent difference measurements
with open("ekg_simulation_output.txt", "a") as f:
    print('\nMean of Percent Difference Measurements:',file=f)
    print('PR-interval:',np.average(pr_error_vec),file=f)
    print('QT-interval:',np.average(qt_error_vec),file=f)
    print('ST-interval:',np.average(st_error_vec),file=f)
    print('QRS-duration:',np.average(qrs_error_vec),file=f)
    print('P-wave duration:',np.average(p_error_vec),file=f)
    print('T-wave duration:',np.average(t_error_vec),file=f)

    print('\nStandard Deviation of Percent Difference Measurements:',file=f)
    print('PR-interval:',np.std(pr_error_vec),file=f)
    print('QT-interval:',np.std(qt_error_vec),file=f)
    print('ST-interval:',np.std(st_error_vec),file=f)
    print('QRS-duration:',np.std(qrs_error_vec),file=f)
    print('P-wave duration:',np.std(p_error_vec),file=f)
    print('T-wave duration:',np.std(t_error_vec),file=f)

