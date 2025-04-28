#!python
from pytmatrix import tmatrix, radar, orientation, refractive, psd
import numpy as np
from pytmatrix import tmatrix_aux

ratio = lambda d: 0.9951 + 0.02510*d - 0.03644*d**2 + 0.005303*d**3 - 0.0002492*d**4
# ratio = lambda d: 1
De = np.logspace(np.log10(0.1),np.log10(8),251)
D = (De[0:-1]+De[1:De.size])/2

wavelengths = ['wl_S', 'wl_C', 'wl_X', 'wl_Ku', 'wl_Ka']
# wavelengths = ['wl_']
# refractive.set_m(wavelengths)

reflh = np.zeros((len(wavelengths),len(D)))
reflv = np.zeros((len(wavelengths),len(D)))
xsecth = np.zeros((len(wavelengths),len(D)))
xsectv = np.zeros((len(wavelengths),len(D)))
rhohva = np.zeros((len(wavelengths),len(D)))
rhohvb = np.zeros((len(wavelengths),len(D)))
kdp = np.zeros((len(wavelengths),len(D)))
# S = np.zeros((len(wavelengths),len(D)))
# M = list(len(wavelengths),len(D))
for wl,i in zip(wavelengths,range(len(wavelengths))):
    print(wl)
    for d,j in zip(D,range(len(D))):
        scatterer = tmatrix.Scatterer(
                                    radius=d/2, 
                                    wavelength=getattr(tmatrix_aux, wl), 
                                    m=refractive.m_w_20C[getattr(tmatrix_aux, wl)], 
                                    axis_ratio=1/ratio(d))
        scatterer.or_pdf = orientation.gaussian_pdf(10.0)
        scatterer.orient = orientation.orient_averaged_adaptive
        scatterer.set_geometry(tmatrix_aux.geom_horiz_back)
        # scatterer.get_S()
        reflh[i,j]=radar.refl(scatterer)
        reflv[i,j]=radar.refl(scatterer,False)
        xsecth[i,j]=radar.radar_xsect(scatterer)
        xsectv[i,j]=radar.radar_xsect(scatterer,False)
        Z = scatterer.get_Z()
        # (Z[2,2] + Z[3,3])**2 + (Z[3,2] - Z[2,3])**2
        rhohva[i,j]=Z[2,2] + Z[3,3]
        rhohvb[i,j]=Z[3,2] - Z[2,3]
        
        scatterer.set_geometry(tmatrix_aux.geom_horiz_forw)
        kdp[i,j] = radar.Kdp(scatterer)
        

np.savetxt('reflh.txt',reflh, delimiter=",")
np.savetxt('reflv.txt',reflv, delimiter=",")
np.savetxt('xsecth.txt',xsecth, delimiter=",")
np.savetxt('xsectv.txt',xsectv, delimiter=",")
np.savetxt('rhohva.txt',rhohva, delimiter=",")
np.savetxt('rhohvb.txt',rhohvb, delimiter=",")
np.savetxt('kdp.txt',kdp, delimiter=",")
