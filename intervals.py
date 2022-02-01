#This script includes functions used to measure the PR-interval, QT-interval, ST-segment, QRS-duration, P-wave duration, and T-wave duration

import numpy as np
from processing import *
from cycles import *



def get_rr_interval(r_peak_xcs):
    #input 1d array of time-coordinates of R-wave peaks
    #output the average RR-interval and standard deviation of this measurement

    rr_int = []

    for i in range(0,len(r_peak_xcs)-1):
        rr_int.append(r_peak_xcs[i+1] - r_peak_xcs[i])

    return np.average(rr_int), np.std(rr_int)



def get_pr_interval(left_p_wave_coords,qrs_onset_coords,r_peak_xcs):
    #input 1d array of P-wave onset time coordinates, 1d array of Q-wave onset time coordinates, and 1d array of R-wave peak time coordinates
    #output PR-interval measurement

    pr_int = []

    for i in range(0,len(r_peak_xcs)-1):
        for j in range(0,len(left_p_wave_coords)):
            for k in range(0,len(qrs_onset_coords)):
                if (left_p_wave_coords[j] > r_peak_xcs[i] and left_p_wave_coords[j] < r_peak_xcs[i+1] and
                qrs_onset_coords[k] > r_peak_xcs[i] and qrs_onset_coords[k] < r_peak_xcs[i+1]):
                    pr_int.append(qrs_onset_coords[k] - left_p_wave_coords[j])

    return np.average(pr_int), np.std(pr_int)



def get_qrs_duration(qrs_onset_coords,right_s_wave_coords,r_peak_xcs):
    #input 1d array of Q-wave onset time coordinates, 1d array of S-wave offset time coordinates, and 1d array of R-wave peak time coordinates
    #output QRS-duration measurement

    qrs_dur = []

    for i in range(0,len(r_peak_xcs)-2):
        for j in range(0,len(qrs_onset_coords)):
            for k in range(0,len(right_s_wave_coords)):
                if (qrs_onset_coords[j] > r_peak_xcs[i] and qrs_onset_coords[j] < r_peak_xcs[i+1]
                and right_s_wave_coords[k] > r_peak_xcs[i+1] and right_s_wave_coords[k] < r_peak_xcs[i+2]):
                    qrs_dur.append(right_s_wave_coords[k] - qrs_onset_coords[j])

    return np.average(qrs_dur), np.std(qrs_dur)



def get_qt_interval(qrs_onset_coords,right_t_wave_coords,r_peak_xcs):
    #input 1d array of Q-wave onset time coordinates, 1d array of T-wave offset time coordinates, and 1d array of R-wave peak time coordinates
    #output QT-interval measurement

    qt_int = []

    for i in range(0,len(r_peak_xcs)-1):
        for j in range(0,len(qrs_onset_coords)):
            for k in range(0,len(right_t_wave_coords)):
                if (qrs_onset_coords[j] > r_peak_xcs[i] and qrs_onset_coords[j] < r_peak_xcs[i+1]
                and right_t_wave_coords[k] > r_peak_xcs[i+1] and right_t_wave_coords[k] < r_peak_xcs[i+2]):
                    qt_int.append(right_t_wave_coords[k] - qrs_onset_coords[j])

    return np.average(qt_int), np.std(qt_int)



def get_st_segment(right_s_wave_coords,left_t_wave_coords,r_peak_xcs):
    #input 1d array of S-wave offset time coordinates, 1d array of T-wave onset time coordinates, and 1d array of R-wave peak time coordinates
    #output ST-segment measurement

    st_int = []
    
    for i in range(0,len(r_peak_xcs)-1):
        for j in range(0,len(right_s_wave_coords)):
            for k in range(0,len(left_t_wave_coords)):
                if (right_s_wave_coords[j] > r_peak_xcs[i] and right_s_wave_coords[j] < r_peak_xcs[i+1] and
                left_t_wave_coords[k] > r_peak_xcs[i] and left_t_wave_coords[k] < r_peak_xcs[i+1]):
                    st_int.append(left_t_wave_coords[k] - right_s_wave_coords[j])

    return np.average(st_int), np.std(st_int)


def get_wave_duration(left_coords,right_coords):
    #input 1d array of left-most points of detected P,T-waves and 1d array of right-most points of detected P,T-waves
    #output the duration of the P or T-wave

    wave_dur = []

    if len(left_coords) == len(right_coords):
        for i in range(0,len(left_coords)):
            wave_dur.append(right_coords[i] - left_coords[i])
        
    return np.average(wave_dur), np.std(wave_dur)


def get_intervals_and_H1wave_idxs(ekg,r_peak_xcs,rr_int_avg,persist,births,xcs,ycs,bps):
    #input ekg signal as nx2 matrix, 1d array as R-wave peak time coordinates, average RR-interval of the signal,
    #1d array of persistent values for 1-dimensional homology features, 1d array of birth radius values for each 1-dimensional homology feature
    #1d array of time-coordinates of representative cycles, 1d array of amplitude-coordinates of representative cycles, and list of boundary points of representative cycles
    #output:12x1 array of means and standard deviations of interval measurements, 4 1d arrays with index locations of P,Q,S, and T-waves in the vector persist

    persist_p_wave_lower = 0.001
    persist_p_wave_upper = 0.2
    persist_t_wave_lower = 0.01
    persist_t_wave_upper = 0.6
    persist_q_wave_lower = 0.007
    persist_q_wave_upper = 0.1
    persist_s_wave_lower = 0.007
    persist_s_wave_upper = 0.1
    p_wave_persist = []
    t_wave_persist = []
    q_wave_persist = []
    s_wave_persist = []

    for i in range(0,len(persist)):
        for j in range(0,len(r_peak_xcs)):
            if (persist[i] > persist_p_wave_lower and persist[i] < persist_p_wave_upper and
            xcs[i] < r_peak_xcs[j]-0.06*rr_int_avg and xcs[i] > r_peak_xcs[j]-0.35*rr_int_avg and
            ycs[i] < 0.15 and ycs[i] > 0.0 and births[i] < 0.03):
                p_wave_persist.append(persist[i])

            if (persist[i] > persist_t_wave_lower and persist[i] < persist_t_wave_upper and
            xcs[i] > r_peak_xcs[j]+0.15*rr_int_avg and xcs[i] < r_peak_xcs[j]+0.5*rr_int_avg and
            ycs[i] < 0.4 and ycs[i] > 0.0 and births[i] < 0.04):
                t_wave_persist.append(persist[i])

            if (persist[i] > persist_q_wave_lower and persist[i] < persist_q_wave_upper and
            xcs[i] < r_peak_xcs[j] and xcs[i] > r_peak_xcs[j]-0.12*rr_int_avg and
            ycs[i] < 0.0 and births[i] < 0.06):
                q_wave_persist.append(persist[i])

            if (persist[i] > persist_s_wave_lower and persist[i] < persist_s_wave_upper and
            xcs[i] > r_peak_xcs[j] and xcs[i] < r_peak_xcs[j]+0.12*rr_int_avg and
            ycs[i] < 0.0 and births[i] < 0.06):
                s_wave_persist.append(persist[i])

    idx_p = []
    idx_t = []
    idx_q = []
    idx_s = []

    for j in range(0,len(persist)):
        for i in range(0,len(p_wave_persist)):
            if persist[j] == p_wave_persist[i]:
                idx_p.append(j)

    for j in range(0,len(persist)):
        for i in range(0,len(t_wave_persist)):
            if persist[j] == t_wave_persist[i]:
                idx_t.append(j)

    for j in range(0,len(persist)):
        for i in range(0,len(q_wave_persist)):
            if persist[j] == q_wave_persist[i]:
                idx_q.append(j)

    for j in range(0,len(persist)):
        for i in range(0,len(s_wave_persist)):
            if persist[j] == s_wave_persist[i]:
                idx_s.append(j)

    left_p_wave_coords = get_left_p_wave_coords(idx_p,bps,r_peak_xcs)
    right_p_wave_coords = get_right_p_wave_coords(idx_p,bps,r_peak_xcs)
    qrs_onset_coords = get_qrs_onset_coords(idx_q,bps,r_peak_xcs,ekg)
    right_s_wave_coords = get_right_s_wave_coords(idx_s,bps,r_peak_xcs)
    left_t_wave_coords = get_left_t_wave_coords(idx_t,bps,r_peak_xcs)
    right_t_wave_coords = get_right_t_wave_coords(idx_t,bps,r_peak_xcs)
    #all_right_s_wave_coords=get_all_right_s_wave_coords(right_s_wave_coords,r_peak_xcs,ekg)

    pr_int_avg,pr_int_sd = get_pr_interval(left_p_wave_coords,qrs_onset_coords,r_peak_xcs)
    qrs_dur_avg,qrs_dur_sd = get_qrs_duration(qrs_onset_coords,right_s_wave_coords,r_peak_xcs)
    qt_int_avg,qt_int_sd = get_qt_interval(qrs_onset_coords,right_t_wave_coords,r_peak_xcs)
    st_int_avg,st_int_sd = get_st_segment(right_s_wave_coords,left_t_wave_coords,r_peak_xcs)
    p_wave_dur_avg,p_wave_dur_sd = get_wave_duration(left_p_wave_coords,right_p_wave_coords)
    t_wave_dur_avg,t_wave_dur_sd = get_wave_duration(left_t_wave_coords,right_t_wave_coords)
    int_tda = [pr_int_avg,pr_int_sd,qt_int_avg,qt_int_sd,st_int_avg,st_int_sd,qrs_dur_avg,qrs_dur_sd,p_wave_dur_avg,p_wave_dur_sd,t_wave_dur_avg,t_wave_dur_sd]
    
    return int_tda, idx_p, idx_q, idx_s, idx_t
