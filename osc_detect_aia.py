from matplotlib import pyplot as plt

import osc_detect 
import numpy as np
import sunpy

def main():
    directory = '/home/ireland/Data/AIA_Data/test2/'
    
    mc = sunpy.make_map(directory).derotate_by_center_of_fov()
    
    return mc
    
if __name__ == '__main__':
    main()