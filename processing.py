#This script includes functions used to process an ECG signal prior to computing its persistent homology

import numpy as np
import statistics as stats



def add_points(data):
    #input a 1d array of length n
    #output a 1d array of length 2n by including the average of each pair of neighboring points as a point in the array
    temp = np.zeros((len(data)))
    temp[len(temp)-1] = data[len(temp)-1]
    data_f = np.zeros((2*len(temp)))

    for i in range(0,len(temp)-1):
        temp[i] = (data[i]+data[i+1])/2
    
    for i in range(0,len(data_f)):
        if i % 2 == 0:
            data_f[i] = data[int(i/2)]
        if i % 2 == 1:
            data_f[i] = temp[int((i-1)/2)]

    return data_f



def adjoin_noise(data, noise):
    #input 1d array data and 1d array noise of same length n
    #returns 1d array of data with noise of length n
    data_f = np.zeros((len(data)))

    for i in range(0,len(data_f)):
        if i % 2 == 0:
            data_f[i] = data[int(i/2)]
        if i % 2 == 1:
            data_f[i] = data[int((i-1)/2)] + noise[int((i-1)/2)]

    return data_f



def reduce_number_points(data):
    #inputs a 1d array and removes every other point
    data_f = []

    for i in range(0,len(data)):
        if i % 2 == 0:
            data_f.append(data[i])

    return(np.array(data_f))
    


def normalize(data):
    #inputs a 1d array, then normalizes it such that its largest element is 1
    max_element = np.max(data)
    data = data / max_element
    return(data)



def add_time_axis(data,sf):
    #input 1d array of voltage measurements of length n and sampling frequency
    #returns n x 2 matrix such that first column is time and second column is voltage
    total_time = len(data) / sf
    ekg = np.zeros(shape=(len(data),2))

    for i in range(0,len(data)):
        for j in range(0,2):
            if j == 0:
                ekg[i][j] = i * total_time / len(data)
            if j == 1:
                ekg[i][j] = data[i]

    return(ekg)



def get_rpeak_xcs(rpeaks,ekg):
    #input 1d array of R-wave peak coordinates and nx2 matrix of EKG signal
    #output 1d array of time coordinates of R-wave peaks 
    r_peak_xc_idx = rpeaks[0]
    r_peak_xcs = []

    for i in range(0,len(ekg)):
        for j in range(0,len(r_peak_xc_idx)):
            if i == r_peak_xc_idx[j]:
                r_peak_xcs.append(ekg[r_peak_xc_idx[j],0])

    r_peak_xcs = np.asarray(r_peak_xcs)

    return r_peak_xcs, r_peak_xc_idx



def add_isoelectric_line(ekg):
    voltage_rounded = np.zeros((len(ekg[:,0])))

    for i in range(0,len(ekg[:,0])):
        voltage_rounded[i] = round(ekg[i,1],2)

    baseline = stats.median(voltage_rounded)

    for i in range(0,len(ekg[:,0])):
        if i % 2 == 0:
            ekg[i,1] = baseline
        if i % 2 == 1:
            ekg[i,1] = ekg[i,1]

    ekg[:,1] = ekg[:,1] - baseline

    return ekg



def add_isoelectric_line2(ekg,r_peak_xc_idx):#this should account for baseline wander, but needs to be modified
    #this function is not currently being used
    #goal of this function is to add an isoelectric line to an ECG signal that accounts for baseline wander
    for j in range(0,len(r_peak_xc_idx)-1):
        voltage_rounded=np.zeros((len(ekg[r_peak_xc_idx[j]:r_peak_xc_idx[j+1],0])))
        voltage_rounded=[]
        for i in range(r_peak_xc_idx[j],r_peak_xc_idx[j+1]):
            voltage_rounded.append(round(ekg[i,1],2))
        baseline=stats.median(voltage_rounded)
        for i in range(r_peak_xc_idx[j],r_peak_xc_idx[j+1]):
            if i%2== 0:
                ekg[i,1] = baseline
        for i in range(r_peak_xc_idx[j],r_peak_xc_idx[j+1]):
            ekg[i,1]=ekg[i,1]-baseline
    return(ekg)



def trim(ekg,r_peak_xcs):
    #input nx2 matrix representing EKG signal and 1d array of R-wave peak time coordinates
    #output mx2 matrix representing EKG signal that starts and ends with an R-wave peak

    for i in range(0,len(ekg[:,0])):
        if ekg[i,0] < r_peak_xcs[0]:
            ekg[i,0] = np.NAN
        if ekg[i,0] > r_peak_xcs[len(r_peak_xcs)-1]:
            ekg[i,0] = np.NAN

    ekg_t_time = [x for x in ekg[:,0] if np.isnan(x) == False]
    ekg_t_volt = np.zeros((len(ekg_t_time)))

    for i in range(0,len(ekg)):
        for j in range(0,len(ekg_t_time)):
            if ekg[i,0] == ekg_t_time[j]:
                ekg_t_volt[j] = ekg[i,1]

    ekg_f = np.zeros((len(ekg_t_time),2))
    ekg_f[:,0] = ekg_t_time
    ekg_f[:,1] = ekg_t_volt
    ekg_f = np.asarray(ekg_f)

    return ekg_f
