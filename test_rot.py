from sunpy.wcs import convert_hg_hcc, convert_hcc_hg
from sunpy.coords import pb0r
from sunpy.coords import rot_hcc

dt = '2008-04-12'
v = pb0r(dt)

x = 200
y = 100
print '============'
print 'Arcsecond location', x, y
h1,h2 = convert_hcc_hg(v["sd"]/60.0, v["b0"], v["l0"], x/3600.0, y/3600.0)
print 'Lat, lon', h1, h2

nx, ny = convert_hg_hcc(v["sd"]/60.0, v["b0"], v["l0"], h1, h2)
print 'Roundtrip arcsec loc.', 3600*nx, 3600*ny

#nx, ny = convert_hg_hcc(v["sd"]/60.0, v["b0"], v["l0"], h1, h2)
#print nx*3600, ny*3600
nnx, nny = convert_hg_hcc(v["sd"]/60.0, v["b0"], v["l0"]   ,    10.931230, 48.717380)
print nnx*3600, ny*3600

newx, newy = rot_hcc(x,y, tstart=dt, tend='2008-04-14')
print newx, newy