from homcloud.interface import optimal_volume
import numpy as np
import homcloud.interface as hc
from numpy.random.mtrand import normal
import pandas as pd
from scipy import signal
from processing import *
from cycles import *
from intervals import *


data=np.asarray(pd.read_csv('/home/hunter/ekg/data/ekg_sim_data.csv',delimiter=','))
intervals=np.asarray(pd.read_csv('/home/hunter/ekg/data/ekg_sim_intervals.csv',delimiter=','))#order: pr,qt,st,qrs,p,t
noise=np.random.normal(0,0.008,len(data[0,:]))
sf=500
x=np.linspace(0,4.0,2*len(data[0,:]))

pr_sim=[]
qt_sim=[]
st_sim=[]
qrs_sim=[]
p_sim=[]
t_sim=[]
int_tda=np.zeros((len(data[:,0]),12))
n_samples=2

pp = PdfPages('ekg_simulation_plots.pdf')
#for i in range(0,len(data[:,0])-1):
for i in range(0,n_samples):
    print(i)
    count=i+1
    data_temp=data[i,:]
    data_temp=add_points(data_temp)
    #data_temp=adjoin_noise(data_temp,noise)
    ekg=np.zeros((len(data_temp),2))
    ekg[:,0]=x[:]
    ekg[:,1]=data_temp[:]
    ekg[:,1]=normalize(ekg[:,1]) #first normalize signal to find R-wave peaks
    rpeaks=signal.find_peaks(ekg[:,1],height=0.5,distance=10)
    r_peak_xcs,r_peak_xc_idx=get_rpeak_xcs(rpeaks,ekg)
    ekg=add_isoelectric_line(ekg,r_peak_xc_idx)
    #ekg[:,1]=adjoin_noise(ekg[:,1],noise)
    #ekg[:,1]=normalize(ekg[:,1]) # now normalize data so that persistent thresholds can be set
    #rpeaks=signal.find_peaks(ekg[:,1],height=0.6,distance=10)
    #r_peak_xcs,r_peak_xc_idx=get_rpeak_xcs(rpeaks,ekg)
    rr_int_avg,rr_int_sd=get_rr_interval(r_peak_xcs)
    ekg=trim(ekg,r_peak_xcs)


    #compute persistent homology
    output=hc.PDList.from_alpha_filtration(ekg,no_squared=True,save_boundary_map=True,save_phtrees=True,save_to="pointcloud.pdgm")
    pd1=hc.PDList("pointcloud.pdgm").dth_diagram(1)

    persist=np.asarray(pd1.deaths-pd1.births)
    births=np.asarray(pd1.births)
    deaths=np.asarray(pd1.deaths)
    #plt.scatter(pd1.births,pd1.deaths)
    #plt.xlabel('Birth Time')
    #plt.ylabel('Death Time')
    #plt.show()


    #cycle_xcs,cycle_ycs,bps=get_vol_opt_cycle_centroid_coords(persist,pd1)
    cycle_xcs,cycle_ycs,bps=get_card_opt_cycle_centroid_coords(persist,pd1)
    #cycle_xcs,cycle_ycs,bps=get_stab_vol_cycle_centroid_coords(persist,pd1)
    #cycle_xcs,cycle_ycs,bps=get_stab_subvol_cycle_centroid_coords(persist,pd1)
    int_tda[i,:],idx_p,idx_q,idx_s,idx_t=get_intervals_and_H1wave_idxs(ekg,r_peak_xcs,rr_int_avg,persist,births,cycle_xcs,cycle_ycs,bps)
    #get_intervals_draw_cycles(ekg,r_peak_xcs,rr_int_avg,persist,births,stab_vol_cycle_xcs,stab_vol_cycle_ycs,bps)
    #get_intervals_draw_cycles(ekg,r_peak_xcs,rr_int_avg,persist,births,stab_subvol_cycle_xcs,stab_subvol_cycle_ycs,bps)
    #print(intervals[temp,:])
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
    fig=draw_cycles(ekg,idx_p,idx_q,idx_s,idx_t,r_peak_xcs,cycle_xcs,bps,count)
    pp.savefig(fig)
pp.close()


pr_tda=[]
qt_tda=[]
st_tda=[]
qrs_tda=[]
p_tda=[]
t_tda=[]
#for i in range(0,len(int[:,0])-1):
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



with open("ekg_simulation_output.txt", "a") as f:
    print('\nAverage Percent Difference in Measurements:',file=f)
    print('PR-interval:',np.average(pr_error_vec),file=f)
    print('QT-interval:',np.average(qt_error_vec),file=f)
    print('ST-interval:',np.average(st_error_vec),file=f)
    print('QRS-duration:',np.average(qrs_error_vec),file=f)
    print('P-wave duration:',np.average(p_error_vec),file=f)
    print('T-wave duration:',np.average(t_error_vec),file=f)




#plt.scatter(opt_vol_xc[len(persist)-8],opt_vol_yc[len(persist)-8],c='orange',marker='.')
#plt.scatter(stable_volume_boundary_points[:,0],stable_volume_boundary_points[:,1],c='orange',marker='p')
#plt.scatter(stable_volume_boundary_points[:,0],stable_subvolume_boundary_points[:,1],c='orange',marker='p')

