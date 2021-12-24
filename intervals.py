from homcloud.interface import optimal_volume
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from processing import *
from cycles import *



def get_rr_interval(r_peak_xcs):
    #input a numpy array of time-coordinates of R-wave peaks
    #output the average RR-interval and standard deviation of this measurement
    rr_int=[]
    for i in range(0,len(r_peak_xcs)-1):
        rr_int.append(r_peak_xcs[i+1]-r_peak_xcs[i])
    return np.average(rr_int),np.std(rr_int)



def get_pr_interval(left_p_wave_coords,qrs_onset_coords,r_peak_xcs):
    pr_int=[]
    for i in range(0,len(r_peak_xcs)-1):
        for j in range(0,len(left_p_wave_coords)):
            for k in range(0,len(qrs_onset_coords)):
                if (left_p_wave_coords[j]>r_peak_xcs[i] and left_p_wave_coords[j]<r_peak_xcs[i+1] and
                qrs_onset_coords[k]>r_peak_xcs[i] and qrs_onset_coords[k]<r_peak_xcs[i+1]):
                    pr_int.append(qrs_onset_coords[k]-left_p_wave_coords[j])
    return np.average(pr_int),np.std(pr_int)



def get_qrs_duration(qrs_onset_coords,right_s_wave_coords,r_peak_xcs):
    qrs_dur=[]
    for i in range(0,len(r_peak_xcs)-2):
        for j in range(0,len(qrs_onset_coords)):
            for k in range(0,len(right_s_wave_coords)):
                if (qrs_onset_coords[j]>r_peak_xcs[i] and qrs_onset_coords[j]<r_peak_xcs[i+1]
                and right_s_wave_coords[k]>r_peak_xcs[i+1] and right_s_wave_coords[k]<r_peak_xcs[i+2]):
                    qrs_dur.append(right_s_wave_coords[k]-qrs_onset_coords[j])
    return np.average(qrs_dur),np.std(qrs_dur)



def get_qt_interval(qrs_onset_coords,right_t_wave_coords,r_peak_xcs):
    qt_int=[]
    for i in range(0,len(r_peak_xcs)-1):
        for j in range(0,len(qrs_onset_coords)):
            for k in range(0,len(right_t_wave_coords)):
                if (qrs_onset_coords[j]>r_peak_xcs[i] and qrs_onset_coords[j]<r_peak_xcs[i+1]
                and right_t_wave_coords[k]>r_peak_xcs[i+1] and right_t_wave_coords[k]<r_peak_xcs[i+2]):
                    qt_int.append(right_t_wave_coords[k]-qrs_onset_coords[j])
    return np.average(qt_int),np.std(qt_int)



def get_st_interval(right_s_wave_coords,left_t_wave_coords,r_peak_xcs):
    st_int=[]
    for i in range(0,len(r_peak_xcs)-1):
        for j in range(0,len(right_s_wave_coords)):
            for k in range(0,len(left_t_wave_coords)):
                if (right_s_wave_coords[j]>r_peak_xcs[i] and right_s_wave_coords[j]<r_peak_xcs[i+1] and
                left_t_wave_coords[k]>r_peak_xcs[i] and left_t_wave_coords[k]<r_peak_xcs[i+1]):
                    st_int.append(left_t_wave_coords[k]-right_s_wave_coords[j])
    return np.average(st_int),np.std(st_int)


def get_wave_duration(left_coords,right_coords):
    wave_dur=[]
    if len(left_coords)==len(right_coords):
        for i in range(0,len(left_coords)):
            wave_dur.append(right_coords[i]-left_coords[i])
        
    return np.average(wave_dur),np.std(wave_dur)


def get_intervals_and_H1wave_idxs(ekg,r_peak_xcs,rr_int_avg,persist,births,cycle_xcs,cycle_ycs,bps):
    persist_tol_p_wave_lower=0.002
    persist_tol_p_wave_upper=0.2
    persist_tol_t_wave_lower=0.01
    persist_tol_t_wave_upper=0.6
    persist_tol_q_wave_lower=0.001
    persist_tol_q_wave_upper=0.08
    persist_tol_s_wave_lower=0.001
    persist_tol_s_wave_upper=0.08
    p_wave_persist=[]
    t_wave_persist=[]
    q_wave_persist=[]
    s_wave_persist=[]
    for i in range(0,len(persist)):
        for j in range(0,len(r_peak_xcs)):
            if (persist[i]>persist_tol_p_wave_lower and persist[i]<persist_tol_p_wave_upper and
            cycle_xcs[i]<r_peak_xcs[j]-0.08*rr_int_avg and cycle_xcs[i]>r_peak_xcs[j]-0.35*rr_int_avg and
            cycle_ycs[i]<0.15 and cycle_ycs[i]>0.0 and births[i]<0.03):
                p_wave_persist.append(persist[i])

            if (persist[i]>persist_tol_t_wave_lower and persist[i]<persist_tol_t_wave_upper and
            cycle_xcs[i]>r_peak_xcs[j]+0.15*rr_int_avg and cycle_xcs[i]<r_peak_xcs[j]+0.6*rr_int_avg and
            cycle_ycs[i]<0.4 and cycle_ycs[i]>0.0 and births[i]<0.08):
                t_wave_persist.append(persist[i])

            if (persist[i]>persist_tol_q_wave_lower and persist[i]<persist_tol_q_wave_upper and
            cycle_xcs[i]<r_peak_xcs[j] and cycle_xcs[i]>r_peak_xcs[j]-0.12*rr_int_avg and
            cycle_ycs[i]<0.0 and births[i]<0.07):
                q_wave_persist.append(persist[i])

            if (persist[i]>persist_tol_s_wave_lower and persist[i]<persist_tol_s_wave_upper and
            cycle_xcs[i]>r_peak_xcs[j] and cycle_xcs[i]<r_peak_xcs[j]+0.12*rr_int_avg and
            cycle_ycs[i]<0.0 and births[i]<0.07):
                s_wave_persist.append(persist[i])

    idx_p=[]
    idx_t=[]
    idx_q=[]
    idx_s=[]
    for j in range(0,len(persist)):
        for i in range(0,len(p_wave_persist)):
            if persist[j]==p_wave_persist[i]:
                idx_p.append(j)
    for j in range(0,len(persist)):
        for i in range(0,len(t_wave_persist)):
            if persist[j]==t_wave_persist[i]:
                idx_t.append(j)
    for j in range(0,len(persist)):
        for i in range(0,len(q_wave_persist)):
            if persist[j]==q_wave_persist[i]:
                idx_q.append(j)
    for j in range(0,len(persist)):
        for i in range(0,len(s_wave_persist)):
            if persist[j]==s_wave_persist[i]:
                idx_s.append(j)

    left_p_wave_coords=get_left_p_wave_coords(idx_p,bps,r_peak_xcs)
    right_p_wave_coords=get_right_p_wave_coords(idx_p,bps,r_peak_xcs)
    qrs_onset_coords=get_qrs_onset_coords(idx_q,bps,r_peak_xcs,ekg)
    right_s_wave_coords=get_right_s_wave_coords(idx_s,bps,r_peak_xcs)
    left_t_wave_coords=get_left_t_wave_coords(idx_t,bps,r_peak_xcs)
    right_t_wave_coords=get_right_t_wave_coords(idx_t,bps,r_peak_xcs)
    all_right_s_wave_coords=get_all_right_s_wave_coords(right_s_wave_coords,r_peak_xcs,ekg)
    #print('\nThe following are the endpoints coming only from detected H1 features:')
    #print('\nThe left-hand endpoints of the P-waves are:\n',left_p_wave_coords)
    #print('\nThe right-hand endpoints of the P-waves are:\n',right_p_wave_coords)
    #print('\nThe right-hand endpoints of the S-waves are:\n',right_s_wave_coords)
    #print('\nThe QRS onset coordinates are:\n',qrs_onset_coords)
    #print('\nThe left-hand endpoints of the T-waves are:\n',left_t_wave_coords)
    #print('\nThe right-hand endpoints of the T-waves are:\n',right_t_wave_coords)
    #print('\n\nAll right-hand S-wave coordinates are:',all_right_s_wave_coords)

    pr_int_avg,pr_int_sd=get_pr_interval(left_p_wave_coords,qrs_onset_coords,r_peak_xcs)
    qrs_dur_avg,qrs_dur_sd=get_qrs_duration(qrs_onset_coords,right_s_wave_coords,r_peak_xcs)
    qt_int_avg,qt_int_sd=get_qt_interval(qrs_onset_coords,right_t_wave_coords,r_peak_xcs)
    st_int_avg,st_int_sd=get_st_interval(right_s_wave_coords,left_t_wave_coords,r_peak_xcs)
    p_wave_dur_avg,p_wave_dur_sd=get_wave_duration(left_p_wave_coords,right_p_wave_coords)
    t_wave_dur_avg,t_wave_dur_sd=get_wave_duration(left_t_wave_coords,right_t_wave_coords)
    int_tda=[pr_int_avg,pr_int_sd,qt_int_avg,qt_int_sd,st_int_avg,st_int_sd,qrs_dur_avg,qrs_dur_sd,p_wave_dur_avg,p_wave_dur_sd,t_wave_dur_avg,t_wave_dur_sd]
    
    return int_tda,idx_p,idx_q,idx_s,idx_t
    #print('\nThe average RR-interval is ',rr_int_avg)
    #print('The standard deviation of our RR-interval measurements is ',rr_int_sd)
    '''
    print('\nThe average PR-interval is ',pr_int_avg)
    print('The standard deviation of our PR-interval measurements is ',pr_int_sd)
    print('\nThe average QT-interval is ',qt_int_avg)
    #print('The corrected QT-interval is ',qt_int_avg/np.sqrt(rr_int_avg))
    print('The standard deviation of our QT-interval measurements is ',qt_int_sd)
    print('\nThe average ST-interval is ',st_int_avg)
    print('The standard deviation of our ST-interval measurements is ',st_int_sd)
    print('\nThe average QRS-duration is ',qrs_dur_avg)
    print('The standard deviation of our QRS-duration measurements is ',qrs_dur_sd)
    print('\nThe average P-wave duration is ',p_wave_dur_avg)
    print('The standard deviation of our P-wave duration measurements is ',p_wave_dur_sd)
    print('\nThe average T-wave duration is ',t_wave_dur_avg)
    print('The standard deviation of our T-wave duration measurements is ',t_wave_dur_sd)
    '''
'''
    with PdfPages('multipage_pdf.pdf') as pdf:
        plt.scatter(ekg[:,0],ekg[:,1],c='blue',marker='.')
        for j in range(0,len(r_peak_xcs)-1):
            temp_xcs=[]
            idx_vec=[]
            for i in range(0,len(idx_p)):
                if (cycle_xcs[idx_p[i]]>r_peak_xcs[j] and cycle_xcs[idx_p[i]]<r_peak_xcs[j+1]):
                    temp_xcs.append(cycle_xcs[idx_p[i]])
                    idx_vec.append(i)
            if len(temp_xcs)>1:
                for i in range(0,len(temp_xcs)):
                    if i==temp_xcs.index(max(temp_xcs)):
                        plt.scatter(bps[idx_p[idx_vec[i]]][:,0],bps[idx_p[idx_vec[i]]][:,1],c='orange',marker='p')
            if len(temp_xcs)==1:
                for i in range(0,len(idx_p)):
                    if i==idx_vec[0]:
                        plt.scatter(bps[idx_p[i]][:,0],bps[idx_p[i]][:,1],c='orange',marker='p')

        for j in range(0,len(r_peak_xcs)-1):
            temp_xcs=[]
            idx_vec=[]
            for i in range(0,len(idx_q)):
                if (cycle_xcs[idx_q[i]]>r_peak_xcs[j] and cycle_xcs[idx_q[i]]<r_peak_xcs[j+1]):
                    temp_xcs.append(cycle_xcs[idx_q[i]])
                    idx_vec.append(i)
            if len(temp_xcs)>1:
                for i in range(0,len(temp_xcs)):
                    if i==temp_xcs.index(max(temp_xcs)):
                        plt.scatter(bps[idx_q[idx_vec[i]]][:,0],bps[idx_q[idx_vec[i]]][:,1],c='green',marker='p')
            if len(temp_xcs)==1:
                for i in range(0,len(idx_q)):
                    if i==idx_vec[0]:
                        plt.scatter(bps[idx_q[i]][:,0],bps[idx_q[i]][:,1],c='green',marker='p')
        
        for j in range(0,len(r_peak_xcs)-1):
            temp_xcs=[]
            idx_vec=[]
            for i in range(0,len(idx_s)):
                if (cycle_xcs[idx_s[i]]>r_peak_xcs[j] and cycle_xcs[idx_s[i]]<r_peak_xcs[j+1]):
                    temp_xcs.append(cycle_xcs[idx_s[i]])
                    idx_vec.append(i)
            if len(temp_xcs)>1:
                for i in range(0,len(temp_xcs)):
                    if i==temp_xcs.index(min(temp_xcs)):
                        plt.scatter(bps[idx_s[idx_vec[i]]][:,0],bps[idx_s[idx_vec[i]]][:,1],c='purple',marker='p')
            if len(temp_xcs)==1:
                for i in range(0,len(idx_s)):
                    if i==idx_vec[0]:
                        plt.scatter(bps[idx_s[i]][:,0],bps[idx_s[i]][:,1],c='purple',marker='p')

        #for i in range(0,len(idx_p)): uncomment this line and the one below to color all detected p-wave cycles
        #    plt.scatter(bps[idx_p[i]][:,0],bps[idx_p[i]][:,1],c='orange',marker='p')
        for i in range(0,len(idx_t)):
            plt.scatter(bps[idx_t[i]][:,0],bps[idx_t[i]][:,1],c='red',marker='p')
        for i in range(0,len(idx_q)):
            plt.scatter(bps[idx_q[i]][:,0],bps[idx_q[i]][:,1],c='green',marker='p')
        for i in range(0,len(idx_s)): 
            plt.scatter(bps[idx_s[i]][:,0],bps[idx_s[i]][:,1],c='purple',marker='p')
        #for i in range(0,len(persist)):
        #    plt.scatter(opt_vol_bps[i][:,0],opt_vol_bps[i][:,1],c='red',marker='p')
        plt.title('EKG Simulation')
        pdf.savefig()
        plt.close()'''
    #plt.show()

    


