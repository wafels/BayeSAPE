from matplotlib import pyplot as plt

import osc_detect
import numpy as np

def main():
    nt = 200
    start = 104
    
    dc = np.random.normal(size=(100,101,nt))
    index = 12.0*np.arange(start, start+nt+1)
    for i in range(0,nt):
        dc[10,10,i] = dc[10,10,i] + np.sin(2*index[i]*np.pi/180.0)
    dcmean = np.mean(dc,axis=-1)
    for i in range(0,nt):
        dc[...,i] = dc[...,i] - dcmean[...]
    plt.plot(dc[10,10,:])
    
    return osc_detect.BTS(dc,index=index).get_pdf()
    
if __name__ == '__main__':
    main()