import numpy as np
import statistics as stats


def add_points(data_temp):
    temp=np.zeros((len(data_temp)))
    temp[len(temp)-1]=data_temp[len(temp)-1]
    data=np.zeros((2*len(temp)))

    for i in range(0,len(temp)-1):
        temp[i]=(data_temp[i]+data_temp[i+1])/2
    
    for i in range(0,len(data)):
        if i%2==0:
            data[i]=data_temp[int(i/2)]
        if i %2==1:
            data[i]=temp[int((i-1)/2)]
    return data



    data=np.zeros((2*len(data_temp)))
    for i in range(0,len(data)):
        if i%2==0:
            data[i]=data_temp[int(i/2)]
        if i%2==1:
            data[i]=data_temp[2*i]
def adjoin_noise(data,noise):
    #input 1d array data and 1d array noise of same length
    #returns 1d array of data with noise
    dataf=np.zeros((len(data)))
    for i in range(0,len(dataf)):
        if i%2==0:
            dataf[i]=data[int(i/2)]
        if i%2==1:
            dataf[i]=data[int((i-1)/2)]+noise[int((i-1)/2)]
    return dataf


def reduce_number_points(data):
    #inputs a 1d array and removes every other point
    data_temp=[]
    for i in range(0,len(data)):
        if i%2==0:
            data_temp.append(data[i])
    return(np.array(data_temp))
    

def normalize(data):
    #inputs a 1d array, then normalizes it such that its largest element is 1.0
    max_element=np.max(data)
    data=data/max_element
    return(data)


def add_time_axis(data,sf):
    #input 1d array of voltage measurements of length n and sampling frequency
    #returns n x 2 matrix such that first column is time and second column is voltage
    total_time=len(data)/sf
    ekg=np.zeros(shape=(len(data),2))
    for i in range(0,len(data)):
        for j in range(0,2):
            if j ==0:
                ekg[i][j]=i*total_time/len(data)
            if j ==1:
                ekg[i][j]=data[i]
    return(ekg)


def get_rpeak_xcs(rpeaks,ekg):
    r_peak_xc_idx=rpeaks[0]
    r_peak_xcs=[]
    for i in range(0,len(ekg)):
        for j in range(0,len(r_peak_xc_idx)):
            if i==r_peak_xc_idx[j]:
                r_peak_xcs.append(ekg[r_peak_xc_idx[j],0])
    r_peak_xcs=np.asarray(r_peak_xcs)
    return r_peak_xcs,r_peak_xc_idx


def add_isoelectric_line(ekg,r_peak_xc_idx):
    voltage_rounded=np.zeros((len(ekg[:,0])))
    for i in range(0,len(ekg[:,0])):
        voltage_rounded[i]=round(ekg[i,1],2)
    baseline=stats.median(voltage_rounded)
    for i in range(0,len(ekg[:,0])):
        if i%2==0:
            ekg[i,1]=baseline
        if i%2==1:
            ekg[i,1]=ekg[i,1]
    ekg[:,1]=ekg[:,1]-baseline
    return ekg




def add_isoelectric_line2(ekg,r_peak_xc_idx):#this should account for baseline wander, but needs to be modified
    for j in range(0,len(r_peak_xc_idx)-1):
        voltage_rounded=np.zeros((len(ekg[r_peak_xc_idx[j]:r_peak_xc_idx[j+1],0])))
        voltage_rounded=[]
        for i in range(r_peak_xc_idx[j],r_peak_xc_idx[j+1]):
            voltage_rounded.append(round(ekg[i,1],2))
        baseline=stats.median(voltage_rounded)
        for i in range(r_peak_xc_idx[j],r_peak_xc_idx[j+1]):
            if i%2==0:
                ekg[i,1]=baseline
        for i in range(r_peak_xc_idx[j],r_peak_xc_idx[j+1]):
            ekg[i,1]=ekg[i,1]-baseline
    return(ekg)


def trim(ekg,r_peak_xcs):#this trims the data such that it starts and ends with an R-wave
    for i in range(0,len(ekg[:,0])):
        if ekg[i,0]<r_peak_xcs[0]:
            ekg[i,0]=np.NAN
        if ekg[i,0]>r_peak_xcs[len(r_peak_xcs)-1]:
            ekg[i,0]=np.NAN
    ekg_t_time = [x for x in ekg[:,0] if np.isnan(x) == False]
    ekg_t_volt=np.zeros((len(ekg_t_time)))
    for i in range(0,len(ekg)):
        for j in range(0,len(ekg_t_time)):
            if ekg[i,0]==ekg_t_time[j]:
                ekg_t_volt[j]=ekg[i,1]
    ekg_f=np.zeros((len(ekg_t_time),2))
    ekg_f[:,0]=ekg_t_time
    ekg_f[:,1]=ekg_t_volt
    ekg_f=np.asarray(ekg_f)
    return ekg_f


'''
def NLM_1dDarbon(signal,Nvar,P,PatchHW):
    if isinstance(P,int): # scalar has been entered; expand into patch sample index vector
        P = P-1 #Python start index from 0
        Pvec = np.array(range(-P,P+1))
    else:
        Pvec = P # use the vector that has been input
    signal = np.array(signal)
    #debug = [];
    N = len(signal)

    denoisedSig = np.empty(len(signal)) #NaN * ones(size(signal));
    denoisedSig[:] = np.nan
    # to simpify, don't bother denoising edges
    iStart = PatchHW+1
    iEnd = N - PatchHW
    denoisedSig[iStart: iEnd] = 0

    #debug.iStart = iStart;
    #debug.iEnd = iEnd;

    # initialize weight normalization
    Z = np.zeros(len(signal))
    cnt = np.zeros(len(signal))

    # convert lambda value to  'h', denominator, as in original Buades papers
    Npatch = 2 * PatchHW + 1
    h = 2 * Npatch * Nvar**2

    for idx in Pvec: # loop over all possible differences: s - t
        # do summation over p - Eq.3 in Darbon
        k = np.array(range(N))
        kplus = k + idx
        igood = np.where((kplus >=0) & (kplus < N)) # ignore OOB data; we could also handle it
        SSD = np.zeros(len(k))
        SSD[igood] = (signal[k[igood]] - signal[kplus[igood]])**2
        Sdx = np.cumsum(SSD)

        for ii in range(iStart,iEnd): # loop over all points 's'
            distance = Sdx[ii + PatchHW] - Sdx[ii - PatchHW-1] #Eq 4;this is in place of point - by - point MSE
            # but note the - 1; we want to icnlude the point ii - iPatchHW

            w = math.exp(-distance/h) # Eq 2 in Darbon
            t = ii + idx # in the papers, this is not made explicit

            if t>0 and t<N:
                denoisedSig[ii] = denoisedSig[ii] + w * signal[t]
                Z[ii] = Z[ii] + w
                #cnt[ii] = cnt[ii] + 1
                #print('ii',ii)
                #print('t',t)
                #print('w',w)
                #print('denoisedSig[ii]', denoisedSig[ii])
                #print('Z[ii]',Z[ii])
     # loop over shifts

    # now apply normalization
    denoisedSig = denoisedSig/(Z + sys.float_info.epsilon)
    denoisedSig[0: PatchHW+1] =signal[0: PatchHW+1]
    denoisedSig[ - PatchHW: ] =signal[- PatchHW: ]
    #debug.Z = Z;
    return denoisedSig
'''