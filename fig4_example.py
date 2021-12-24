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
sf=500
x=np.linspace(0,2.0,int(len(data[0,:])/2))


int_tda=np.zeros((12))
data_temp=data[5,0:int(len(data[0,:])/2)]
ekg=np.zeros((len(data_temp),2))
ekg[:,0]=x[:]
ekg[:,1]=data_temp[:]
ekg[:,1]=normalize(ekg[:,1]) #first normalize signal to find R-wave peaks
rpeaks=signal.find_peaks(ekg[:,1],height=0.5,distance=10)
r_peak_xcs,r_peak_xc_idx=get_rpeak_xcs(rpeaks,ekg)
ekg=add_isoelectric_line(ekg,r_peak_xc_idx)
#ekg[:,1]=adjoin_noise(ekg[:,1],noise)
ekg[:,1]=normalize(ekg[:,1]) # now normalize data so that persistent thresholds can be set
rpeaks=signal.find_peaks(ekg[:,1],height=0.6,distance=10)
r_peak_xcs,r_peak_xc_idx=get_rpeak_xcs(rpeaks,ekg)
rr_int_avg,rr_int_sd=get_rr_interval(r_peak_xcs)
ekg=trim(ekg,r_peak_xcs)


#compute persistent homology
output=hc.PDList.from_alpha_filtration(ekg,no_squared=True,save_boundary_map=True,save_phtrees=True,save_to="pointcloud.pdgm")
pd1=hc.PDList("pointcloud.pdgm").dth_diagram(1)

persist=np.asarray(pd1.deaths-pd1.births)
births=np.asarray(pd1.births)
deaths=np.asarray(pd1.deaths)

cycle_xcs,cycle_ycs,bps=get_vol_opt_cycle_centroid_coords(persist,pd1)
int_tda,idx_p,idx_q,idx_s,idx_t=get_intervals_and_H1wave_idxs(ekg,r_peak_xcs,rr_int_avg,persist,births,cycle_xcs,cycle_ycs,bps)

pr_sim=intervals[5,0]
qt_sim=intervals[5,1]
st_sim=intervals[5,2]
qrs_sim=intervals[5,3]
p_sim=intervals[5,4]
t_sim=intervals[5,5]

pr_tda=int_tda[0]
qt_tda=int_tda[2]
st_tda=int_tda[4]
qrs_tda=int_tda[6]
p_tda=int_tda[8]
t_tda=int_tda[10]



fig = plt.figure(figsize=[10,7])
ax1 = fig.add_subplot(1, 3, 1)
ax2 = fig.add_subplot(1, 3, 2)
ax3 = fig.add_subplot(1, 3, 3)

p_legend=0
q_legend=0
s_legend=0
t_legend=0
ax1.scatter(ekg[:,0],ekg[:,1],c='blue',marker='.')
for j in range(0,len(r_peak_xcs)-1):
    temp_xcs=[]
    idx_vec=[]
    for i in range(0,len(idx_p)):
        if (cycle_xcs[idx_p[i]]>r_peak_xcs[j] and cycle_xcs[idx_p[i]]<r_peak_xcs[j+1]):
            temp_xcs.append(cycle_xcs[idx_p[i]])
            idx_vec.append(i)
    if len(temp_xcs)>1:
        for i in range(0,len(temp_xcs)):
            if i==temp_xcs.index(max(temp_xcs)) and p_legend<1:
                p_legend=p_legend+1
                ax1.scatter(bps[idx_p[idx_vec[i]]][:,0],bps[idx_p[idx_vec[i]]][:,1],c='orange',marker='p',label='P-wave')
            if i==temp_xcs.index(max(temp_xcs)) and p_legend>=1:
                ax1.scatter(bps[idx_p[idx_vec[i]]][:,0],bps[idx_p[idx_vec[i]]][:,1],c='orange',marker='p')
    if len(temp_xcs)==1:
        for i in range(0,len(idx_p)):
            if i==idx_vec[0] and p_legend<1:
                p_legend=p_legend+1
                ax1.scatter(bps[idx_p[i]][:,0],bps[idx_p[i]][:,1],c='orange',marker='p',label='P-wave')
            if i==idx_vec[0] and p_legend>=1:
                ax1.scatter(bps[idx_p[i]][:,0],bps[idx_p[i]][:,1],c='orange',marker='p')

for j in range(0,len(r_peak_xcs)-1):
    temp_xcs=[]
    idx_vec=[]
    for i in range(0,len(idx_q)):
        if (cycle_xcs[idx_q[i]]>r_peak_xcs[j] and cycle_xcs[idx_q[i]]<r_peak_xcs[j+1]):
            temp_xcs.append(cycle_xcs[idx_q[i]])
            idx_vec.append(i)
    if len(temp_xcs)>1:
        for i in range(0,len(temp_xcs)):
            if i==temp_xcs.index(max(temp_xcs)) and q_legend<1:
                q_legend=q_legend+1
                ax1.scatter(bps[idx_q[idx_vec[i]]][:,0],bps[idx_q[idx_vec[i]]][:,1],c='green',marker='p',label='Q-wave')
            if i==temp_xcs.index(max(temp_xcs)) and q_legend>=1:
                ax1.scatter(bps[idx_q[idx_vec[i]]][:,0],bps[idx_q[idx_vec[i]]][:,1],c='green',marker='p')
    if len(temp_xcs)==1:
        for i in range(0,len(idx_q)):
            if i==idx_vec[0] and q_legend<1:
                q_legend=q_legend+1
                ax1.scatter(bps[idx_q[i]][:,0],bps[idx_q[i]][:,1],c='green',marker='p',label='Q-wave')
            if i==idx_vec[0] and q_legend>=1:
                ax1.scatter(bps[idx_q[i]][:,0],bps[idx_q[i]][:,1],c='green',marker='p')
    
for j in range(0,len(r_peak_xcs)-1):
    temp_xcs=[]
    idx_vec=[]
    for i in range(0,len(idx_s)):
        if (cycle_xcs[idx_s[i]]>r_peak_xcs[j] and cycle_xcs[idx_s[i]]<r_peak_xcs[j+1]):
            temp_xcs.append(cycle_xcs[idx_s[i]])
            idx_vec.append(i)
    if len(temp_xcs)>1:
        for i in range(0,len(temp_xcs)):
            if i==temp_xcs.index(min(temp_xcs)) and s_legend<1:
                s_legend=s_legend+1
                ax1.scatter(bps[idx_s[idx_vec[i]]][:,0],bps[idx_s[idx_vec[i]]][:,1],c='purple',marker='p',label='S-wave')
            if i==temp_xcs.index(min(temp_xcs)) and s_legend>=1:
                ax1.scatter(bps[idx_s[idx_vec[i]]][:,0],bps[idx_s[idx_vec[i]]][:,1],c='purple',marker='p')
    if len(temp_xcs)==1:
        for i in range(0,len(idx_s)):
            if i==idx_vec[0] and s_legend<1:
                s_legend=s_legend+1
                ax1.scatter(bps[idx_s[i]][:,0],bps[idx_s[i]][:,1],c='purple',marker='p',label='S-wave')
            if i==idx_vec[0] and s_legend>=1:
                ax1.scatter(bps[idx_s[i]][:,0],bps[idx_s[i]][:,1],c='purple',marker='p')

for i in range(0,len(idx_t)):
    if i==0:
        ax1.scatter(bps[idx_t[i]][:,0],bps[idx_t[i]][:,1],c='red',marker='p',label='T-wave')
    else:
        ax1.scatter(bps[idx_t[i]][:,0],bps[idx_t[i]][:,1],c='red',marker='p')
ax1.set_title('EKG Simulation with Cycle Representatives')
ax1.set_xlabel('Time')
ax1.set_ylabel('Amplitude')
ax1.legend()


ax2.scatter(pd1.births,pd1.deaths)
ax2.scatter(pd1.births[idx_p],pd1.deaths[idx_p],c='orange',marker='p')
ax2.scatter(pd1.births[idx_q],pd1.deaths[idx_q],c='green',marker='p')
ax2.scatter(pd1.births[idx_s],pd1.deaths[idx_s],c='purple',marker='p')
ax2.scatter(pd1.births[idx_t],pd1.deaths[idx_t],c='red',marker='p')
ax2.set_xlabel('Birth Time')
ax2.set_ylabel('Death Time')


ax3.text(0,0.8,'Simulation-Measured Intervals:')
ax3.text(0,0.77,'PR-interval: {}'.format(float(round(pr_sim,3))))
ax3.text(0,0.74,'QT-interval: {}'.format(float(round(qt_sim,3))))
ax3.text(0,0.71,'ST-segment: {}'.format(float(round(st_sim,3))))
ax3.text(0,0.68,'QRS-duration: {}'.format(float(round(qrs_sim,3))))
ax3.text(0,0.65,'P-wave duration: {}'.format(float(round(p_sim,3))))
ax3.text(0,0.62,'T-wave duration: {}'.format(float(round(t_sim,3))))

ax3.text(0.65,0.8,'TDA-Measured Intervals:')
ax3.text(0.65,0.77,'PR-interval: {}'.format(float(round(pr_tda,3))))
ax3.text(0.65,0.74,'QT-interval: {}'.format(float(round(qt_tda,3))))
ax3.text(0.65,0.71,'ST-segment: {}'.format(float(round(st_tda,3))))
ax3.text(0.65,0.68,'QRS-duration: {}'.format(float(round(qrs_tda,3))))
ax3.text(0.65,0.65,'P-wave duration: {}'.format(float(round(p_tda,3))))
ax3.text(0.65,0.62,'T-wave duration: {}'.format(float(round(t_tda,3))))
ax3.axis('off')
plt.show()
