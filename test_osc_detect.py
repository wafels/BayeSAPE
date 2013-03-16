import sunpy
import osc_detect
import numpy as np

nt = 200
start = 104

dc = np.random.normal(size=(100,101,nt))
index = 12.0*np.arange(start, start+nt+1)

print osc_detect.BTS(dc,index=index).nmap