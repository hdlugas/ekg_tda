#this script identifies P,Q,S,T-waves for real ekg data
import numpy as np
import pandas as pd
import os, random
import homcloud.interface as hc
from numpy.random.mtrand import normal
from matplotlib.backends.backend_pdf import PdfPages
from scipy import signal
from processing import *
from cycles import *
from intervals import *

n_samples=2
sf=500
name_vec=[]

pp = PdfPages('ekg_plots_vol_opt.pdf')
for i in range(0,n_samples):
    print(i)
    count=i+1
    name=random.choice(os.listdir("/home/hunter/ekg/ECGDataDenoised/"))
    name_vec.append(name)
    data=np.asarray(pd.read_csv('/home/hunter/ekg/ECGDataDenoised/'+name,delimiter=','))
    data_temp=data[0:2000,1]
    data_temp=add_points(data_temp)
    x=np.linspace(0,4.0,len(data_temp))
    #data_temp=adjoin_noise(data_temp,noise)
    ekg=np.zeros((len(data_temp),2))
    ekg[:,0]=x[:]
    ekg[:,1]=data_temp[:]
    ekg[:,1]=normalize(ekg[:,1]) #first normalize signal to find R-wave peaks
    rpeaks=signal.find_peaks(ekg[:,1],height=0.5,distance=10)
    r_peak_xcs,r_peak_xc_idx=get_rpeak_xcs(rpeaks,ekg)
    ekg=add_isoelectric_line(ekg)
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

    temp=np.empty((12))
    cycle_xcs,cycle_ycs,bps=get_vol_opt_cycle_centroid_coords(persist,pd1)
    #cycle_xcs,cycle_ycs,bps=get_card_opt_cycle_centroid_coords(persist,pd1)
    #cycle_xcs,cycle_ycs,bps=get_stab_vol_cycle_centroid_coords(persist,pd1)
    #cycle_xcs,cycle_ycs,bps=get_stab_subvol_cycle_centroid_coords(persist,pd1)
    temp,idx_p,idx_q,idx_s,idx_t=get_intervals_and_H1wave_idxs(ekg,r_peak_xcs,rr_int_avg,persist,births,cycle_xcs,cycle_ycs,bps)
    #get_intervals_draw_cycles(ekg,r_peak_xcs,rr_int_avg,persist,births,stab_vol_cycle_xcs,stab_vol_cycle_ycs,bps)
    #get_intervals_draw_cycles(ekg,r_peak_xcs,rr_int_avg,persist,births,stab_subvol_cycle_xcs,stab_subvol_cycle_ycs,bps)

    fig=draw_cycles(ekg,idx_p,idx_q,idx_s,idx_t,r_peak_xcs,cycle_xcs,bps,count)
    pp.savefig(fig)
pp.close()


with open("ecg_file_names_vol_opt.txt", "a") as f:
    for i in range(0,n_samples):
        print('#',i+1,' '+name_vec[i],file=f)