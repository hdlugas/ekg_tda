#This script includes functions related to computing the centroid of boundary points of optimal cycles, 
#computing the onset and offset of the PR-interval, QT-interval, ST-segment, QRS-duration, P-wave duration, and T-wave duration,
#and drawing optimal 1-cycles identified as P,Q,S, and T-waves on EKG signals

import numpy as np
import homcloud.interface as hc
import matplotlib.pyplot as plt



def get_centroid(bps):
    #input list of boundary points
    #output time-coordinate and amplitude-coordinate of centroid of boundary points
    xc = np.average(bps[:,0])
    yc = np.average(bps[:,1])
    return xc, yc



def get_vol_opt_cycle_centroid_coords(persist,pd1):
    #input list of boundary points
    #output time-coordinate and amplitude-coordinate of centroid of boundary points

    bps = []
    xcs = np.zeros((len(persist)))
    ycs = np.zeros((len(persist)))

    for i in range(0,len(persist)):
        pair = hc.Pair(pd1,i)
        opt_vol = pair.optimal_volume()
        bps.append(np.asarray(opt_vol.boundary_points()))
        xcs[i], ycs[i] = get_centroid(bps[i])

    return xcs, ycs, bps


def get_card_opt_cycle_centroid_coords(persist,pd1):
    #input list of boundary points
    #output time-coordinate and amplitude-coordinate of centroid of boundary points

    bps = []
    xcs = np.zeros((len(persist)))
    ycs = np.zeros((len(persist)))

    for i in range(0,len(persist)):
        pair = hc.Pair(pd1,i)
        opt_vol = pair.optimal_1_cycle()
        bps.append(np.asarray(opt_vol.boundary_points()))
        xcs[i], ycs[i] = get_centroid(bps[i])

    return xcs, ycs, bps


def get_stab_vol_cycle_centroid_coords(persist,pd1):
    #input list of boundary points
    #output time-coordinate and amplitude-coordinate of centroid of boundary points

    bps = []
    xcs = np.zeros((len(persist)))
    ycs = np.zeros((len(persist)))

    for i in range(0,len(persist)):
        pair = hc.Pair(pd1,i)
        stab_vol = pair.stable_volume(pow(10,-11))
        bps.append(np.asarray(stab_vol.boundary_points()))
        xcs[i], ycs[i] = get_centroid(bps[i])

    return xcs, ycs, bps



def get_stab_subvol_cycle_centroid_coords(persist,pd1):
    #input list of boundary points
    #output time-coordinate and amplitude-coordinate of centroid of boundary points

    bps = []
    xcs = np.zeros((len(persist)))
    ycs = np.zeros((len(persist)))

    for i in range(0,len(persist)):
        pair = hc.Pair(pd1,i)
        optimal_volume = pair.optimal_volume()
        stab_subvol = optimal_volume.stable_subvolume(pow(10,-11))
        bps.append(np.asarray(stab_subvol.boundary_points()))
        xcs[i], ycs[i] = get_centroid(bps[i])

    return xcs, ycs, bps


def get_left_p_wave_coords(idx_p,opt_vol_bps,r_peak_xcs): 
    #input: 1d array of P-wave indices in the array persist, list of boundary points, and 1d array of R-wave peak time coordinates
    #output: array of the left-hand endpoint of the right-most detected p-wave (if one exists) in each RR-interval

    left_raw = []
    left = []

    for i in range(0,len(idx_p)):
        left_raw.append(np.min(opt_vol_bps[idx_p[i]][:,0]))

    left_raw = sorted(left_raw,key=float)

    for j in range(0,len(r_peak_xcs)-1):
        temp = []
        for i in range(0,len(left_raw)):
            if left_raw[i] > r_peak_xcs[j] and left_raw[i] < r_peak_xcs[j+1]:
                temp.append(left_raw[i])
        if len(temp) > 0:
            left.append(np.max(temp))

    return left    



def get_right_p_wave_coords(idx_p,opt_vol_bps,r_peak_xcs):
    #input: 1d array of P-wave indices in the array persist, list of boundary points, and 1d array of R-wave peak time coordinates
    #output: array of the left-hand endpoint of the right-most detected P-wave (if one exists) in each RR-interval

    right_raw = []
    right = []

    for i in range(0,len(idx_p)):
        right_raw.append(np.max(opt_vol_bps[idx_p[i]][:,0]))

    right_raw = sorted(right_raw,key=float)

    for j in range(0,len(r_peak_xcs)-1):
        temp = []
        for i in range(0,len(right_raw)):
            if right_raw[i] > r_peak_xcs[j] and right_raw[i] < r_peak_xcs[j+1]:
                temp.append(right_raw[i])
        if len(temp) > 0:
            right.append(np.max(temp))

    return right



def get_qrs_onset_coords(idx_q,opt_vol_bps,r_peak_xcs,ekg):
    #input: 1d array of Q-wave indices in the array persist, list of boundary points, and 1d array of R-wave peak time coordinates
    #output: array of the left-hand endpoint of the left-most detected Q-wave (if one exists) in each RR-interval

    qrs_onset_raw = []
    qrs_onset = []

    for i in range(0,len(idx_q)):
        qrs_onset_raw.append(np.min(opt_vol_bps[idx_q[i]][:,0]))

    qrs_onset_raw = sorted(qrs_onset_raw,key=float)

    for j in range(0,len(r_peak_xcs)-1):
        temp = []
        for i in range(0,len(qrs_onset_raw)):
            if qrs_onset_raw[i] > r_peak_xcs[j] and qrs_onset_raw[i] < r_peak_xcs[j+1]:
                temp.append(qrs_onset_raw[i])
        if len(temp) > 0:
            qrs_onset.append(np.max(temp))

    #consider adding qrs onset points when H1 feature of q-wave is not detected
    return qrs_onset     



def get_right_s_wave_coords(idx_s,opt_vol_bps,r_peak_xcs):
    #input: 1d array of S-wave indices in the array persist, list of boundary points, and 1d array of R-wave peak time coordinates
    #output: array of the right-hand endpoint of the right-most detected S-wave (if one exists) in each RR-interval

    right_raw = []
    right = []

    for i in range(0,len(idx_s)):
        right_raw.append(np.max(opt_vol_bps[idx_s[i]][:,0]))

    right_raw = sorted(right_raw,key=float)

    for j in range(0,len(r_peak_xcs)-1):
        temp = []
        for i in range(0,len(right_raw)):
            if right_raw[i] > r_peak_xcs[j] and right_raw[i] < r_peak_xcs[j+1]:
                temp.append(right_raw[i])
        if len(temp) > 0:
            right.append(np.min(temp))

    return right


def get_left_t_wave_coords(idx_t,opt_vol_bps,r_peak_xcs):
    #input: 1d array of T-wave indices in the array persist, list of boundary points, and 1d array of R-wave peak time coordinates
    #output: array of the left-hand endpoint of the left-most detected T-wave (if one exists) in each RR-interval

    left_raw = []
    left = []

    for i in range(0,len(idx_t)):
        left_raw.append(np.min(opt_vol_bps[idx_t[i]][:,0]))

    left_raw = sorted(left_raw,key=float)

    for j in range(0,len(r_peak_xcs)-1):
        temp = []
        for i in range(0,len(left_raw)):
            if left_raw[i] > r_peak_xcs[j] and left_raw[i] < r_peak_xcs[j+1]:
                temp.append(left_raw[i])
        if len(temp) > 0:
            left.append(np.min(temp))

    return left   



def get_right_t_wave_coords(idx_t,opt_vol_bps,r_peak_xcs):
    #input: 1d array of T-wave indices in the array persist, list of boundary points, and 1d array of R-wave peak time coordinates
    #output: array of the right-hand endpoint of the left-most detected T-wave (if one exists) in each RR-interval

    right_raw = []
    right = []

    for i in range(0,len(idx_t)):
        right_raw.append(np.max(opt_vol_bps[idx_t[i]][:,0]))

    right_raw = sorted(right_raw,key=float)

    for j in range(0,len(r_peak_xcs)-1):
        temp = []
        for i in range(0,len(right_raw)):
            if right_raw[i] > r_peak_xcs[j] and right_raw[i] < r_peak_xcs[j+1]:
                temp.append(right_raw[i])
        if len(temp) > 0:
            right.append(np.min(temp))

    return right 



def get_all_right_s_wave_coords(right_s_wave_coords,r_peak_xcs,ekg):
#finish this function to make it return the first x-coordinate of the first datum from ekg above the baseline after the dip from the r-wave
    all_right_s_wave_coords = np.zeros((len(r_peak_xcs)-1))
    for i in range(0,len(r_peak_xcs)-1):
        for j in range(0,len(right_s_wave_coords)):
            if right_s_wave_coords[j]>r_peak_xcs[i] and right_s_wave_coords[j]<r_peak_xcs[i+1]:
                all_right_s_wave_coords[i]=right_s_wave_coords[j]
            else:
                print(' ')
    return all_right_s_wave_coords


def draw_cycles(ekg,idx_p,idx_q,idx_s,idx_t,r_peak_xcs,cycle_xcs,bps,count):
    #input nx2 matrix as ekg signal, 1d arrays with elements being the index location of P,Q,S, and T-waves in the 1d array persist,
    #1d array of R-wave peak time coordinates, 1d array of time coordinates of cycle representatives, list of boundary points of cycles, and the iteration number
    #output a figure of EKG signal with representative cycles identified as P,Q,S, and T-waves drawn

    #note that in a given RR-interval, this function draws the P,S-wave with the right-most centroid and the Q-wave with the left-most centroid. This is 
    #consistent with the 1-cycles used to measure intervals of interest

    p_legend = 0
    q_legend = 0
    s_legend = 0
    t_legend = 0

    fig = plt.figure()
    plt.scatter(ekg[:,0],ekg[:,1],c='blue',marker='.')


    for j in range(0,len(r_peak_xcs)-1):
        temp_xcs = []
        idx_vec = []

        for i in range(0,len(idx_p)):
            if (cycle_xcs[idx_p[i]] > r_peak_xcs[j] and cycle_xcs[idx_p[i]] < r_peak_xcs[j+1]):
                temp_xcs.append(cycle_xcs[idx_p[i]])
                idx_vec.append(i)

        if len(temp_xcs) > 1:
            for i in range(0,len(temp_xcs)):
                if i == temp_xcs.index(max(temp_xcs)) and p_legend < 1:
                    p_legend = p_legend + 1
                    plt.scatter(bps[idx_p[idx_vec[i]]][:,0],bps[idx_p[idx_vec[i]]][:,1],c='orange',marker='p',label='P-wave')
                if i == temp_xcs.index(max(temp_xcs)) and p_legend >= 1:
                    plt.scatter(bps[idx_p[idx_vec[i]]][:,0],bps[idx_p[idx_vec[i]]][:,1],c='orange',marker='p')

        if len(temp_xcs) == 1:
            for i in range(0,len(idx_p)):
                if i == idx_vec[0] and p_legend < 1:
                    p_legend = p_legend + 1
                    plt.scatter(bps[idx_p[i]][:,0],bps[idx_p[i]][:,1],c='orange',marker='p',label='P-wave')
                if i == idx_vec[0] and p_legend >= 1:
                    plt.scatter(bps[idx_p[i]][:,0],bps[idx_p[i]][:,1],c='orange',marker='p')


    for j in range(0,len(r_peak_xcs)-1):
        temp_xcs = []
        idx_vec = []

        for i in range(0,len(idx_q)):
            if (cycle_xcs[idx_q[i]] > r_peak_xcs[j] and cycle_xcs[idx_q[i]] < r_peak_xcs[j+1]):
                temp_xcs.append(cycle_xcs[idx_q[i]])
                idx_vec.append(i)

        if len(temp_xcs) > 1:
            for i in range(0,len(temp_xcs)):
                if i == temp_xcs.index(max(temp_xcs)) and q_legend < 1:
                    q_legend = q_legend + 1
                    plt.scatter(bps[idx_q[idx_vec[i]]][:,0],bps[idx_q[idx_vec[i]]][:,1],c='green',marker='p',label='Q-wave')
                if i == temp_xcs.index(max(temp_xcs)) and q_legend >= 1:
                    plt.scatter(bps[idx_q[idx_vec[i]]][:,0],bps[idx_q[idx_vec[i]]][:,1],c='green',marker='p')

        if len(temp_xcs) == 1:
            for i in range(0,len(idx_q)):
                if i == idx_vec[0] and q_legend < 1:
                    q_legend = q_legend + 1
                    plt.scatter(bps[idx_q[i]][:,0],bps[idx_q[i]][:,1],c='green',marker='p',label='Q-wave')
                if i == idx_vec[0] and q_legend >= 1:
                    plt.scatter(bps[idx_q[i]][:,0],bps[idx_q[i]][:,1],c='green',marker='p')
        

    for j in range(0,len(r_peak_xcs)-1):
        temp_xcs = []
        idx_vec = []

        for i in range(0,len(idx_s)):
            if (cycle_xcs[idx_s[i]] > r_peak_xcs[j] and cycle_xcs[idx_s[i]] < r_peak_xcs[j+1]):
                temp_xcs.append(cycle_xcs[idx_s[i]])
                idx_vec.append(i)

        if len(temp_xcs) > 1:
            for i in range(0,len(temp_xcs)):
                if i == temp_xcs.index(min(temp_xcs)) and s_legend < 1:
                    s_legend = s_legend + 1
                    plt.scatter(bps[idx_s[idx_vec[i]]][:,0],bps[idx_s[idx_vec[i]]][:,1],c='purple',marker='p',label='S-wave')
                if i == temp_xcs.index(min(temp_xcs)) and s_legend >= 1:
                    plt.scatter(bps[idx_s[idx_vec[i]]][:,0],bps[idx_s[idx_vec[i]]][:,1],c='purple',marker='p')

        if len(temp_xcs) == 1:
            for i in range(0,len(idx_s)):
                if i == idx_vec[0] and s_legend < 1:
                    s_legend = s_legend + 1
                    plt.scatter(bps[idx_s[i]][:,0],bps[idx_s[i]][:,1],c='purple',marker='p',label='S-wave')
                if i == idx_vec[0] and s_legend >= 1:
                    plt.scatter(bps[idx_s[i]][:,0],bps[idx_s[i]][:,1],c='purple',marker='p')


    #for i in range(0,len(idx_p)): uncomment this line and the one below to color all detected p-wave cycles
    #    plt.scatter(bps[idx_p[i]][:,0],bps[idx_p[i]][:,1],c='orange',marker='p')
    for i in range(0,len(idx_t)):
        if i==0:
            plt.scatter(bps[idx_t[i]][:,0],bps[idx_t[i]][:,1],c='red',marker='p',label='T-wave')
        else:
            plt.scatter(bps[idx_t[i]][:,0],bps[idx_t[i]][:,1],c='red',marker='p')
    #uncomment the 4 lines below to draw all detected Q and S waves
    #for i in range(0,len(idx_q)):
    #    plt.scatter(bps[idx_q[i]][:,0],bps[idx_q[i]][:,1],c='green',marker='p')
    #for i in range(0,len(idx_s)): 
    #    plt.scatter(bps[idx_s[i]][:,0],bps[idx_s[i]][:,1],c='purple',marker='p')
    #uncomment the 2 lines below to draw all detected P,Q,S, and T-waves
    #for i in range(0,len(persist)):
    #    plt.scatter(opt_vol_bps[i][:,0],opt_vol_bps[i][:,1],c='red',marker='p')

    #plt.title('EKG Simulation #'+str(count)+' Cycle Representatives')
    plt.title('ECG Data with Cycle Representatives')
    plt.xlabel('Time axis')
    plt.ylabel('Amplitude')
    plt.legend()
    #plt.show()

    return fig
