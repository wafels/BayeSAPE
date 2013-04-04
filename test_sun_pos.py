from sunpy import coords

def main():
    x = 0.0
    y = 0.0
    nx, ny = coords.rot_hcc(x, y, tstart='2012/08/18', tend='2012/08/19')
    print nx, ny

if __name__ == '__main__':
    main() 