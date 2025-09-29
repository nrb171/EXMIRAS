
# `radar.m`
This document describes the `radar` (`ra`) class. This class is a toolkit for simulating radar observables from a known drop size distribution. The solver uses LUTs from T-matrix simulations to retrieve the scattering matrix elements for a given radar band and drop size. The `radar` class is used internally by the `exmiras` and `dsdAssimilation` class to simulate the radar variables. You can also use this class independently to simulate radar variables.

To initialize a `radar` object, use:
```matlab
ra = radar(bandName);
```

## `radar` methods
### `N`, the discrete drop size distribution
Please note, for all of the methods requiring `N`, `N` should be a drop size distribution in units of mm^-1 m^-3, defined on the diameter grid `ra.D` (mm), with bin edges defined by `ra.De` (mm). `N` will be an array of length 250. You can create DSDs using the `dsdAssimilation` class or by using your own methods.

To convert between an arbitrary DSD and the format required by the `radar` class, you must calculate the integral of the arbitrary DSD from the edges of the diameter bins (`ra.De`), and divide by the bin width (`ra.Dw`) to convert back to units of mm^-1 m^-3:
$$N_i = \frac{1}{D_{w,i}}\int_{D_{e, i-1}}^{D_{e, i}}n_a(D)dD.$$
Where $N_i$ is the $i$ th element of $N$, and $n_a(D)$ is the arbitrary DSD as a function of diameter.

For your convenience, the `radar` class automatically calculates `ra.D`, `ra.De`, and `ra.Dw` when the object is initialized:
```matlab
ra.De = logspace(log10(0.1),log10(8),ra.nBins+1); % logarithmically spaced bin edges from 0.1 to 8 mm
ra.D =  (ra.De(1:end-1) + ra.De(2:end))/2; % bin centers
ra.Dw = diff(ra.De); % bin widths
```
The `ra.D`, `ra.Dw`, and `ra.De` properties cannot be changed, since the scattering matrix LUTs are defined on this diameter grid. 

### T-matrix LUTs
The `radar` class uses T-matrix scattering matrix element lookup LUTs. The code used to generate the LUTs using the [`pytmatrix`](https://github.com/jleinonen/pytmatrix) high-level interface is available at [tmatrixInterface.py](../scripts/helpers/tmatrixInterface.py). The LUTs were calculated for S, C, X, W, and Ka bands (wavelengths of 111, 53.5, 33.3, 22.0, and 8.43 mm respectively). While `radar` can use any of these wavelengths, only S, C, and X bands have been tested in `exmiras` and `dsdAssimilation` classes, so there may be some issues when using the smallest wavelengths. The droplets T-matrix method was run over 50 droplets with an average canting angle of 0 $^\circ$ and standard deviation of 10 $^\circ$ the scattering matrices were averaged and saved. The beam is assumed to be horizontal [^1].

[^1]: This is a simple implementation of the T-matrix, there is no consideration of horizontal wind shear or turbulence effects on the canting angle distribution, attenuation effects, refractive effects, or different elevation angles. For a more robust implementation of these methodologies, please see [CR-SIM](https://github.com/marikooue/CR-SIM).


 The LUTs are stored in the `EXMIRAS/data/LUTs/` directory.


### Radar variable calculation methods
All of these methods will use an input `N` (mm^-1 m^-3; 1 x 250 array) and return the radar variable in the appropriate units.
- `[Zhh] = ra.calcZhh(N)`: Given a drop size distribution `N` (mm^-1 m^-3), estimate the equivalent horizontal radar reflectivity `Zhh` in units of dBZ. 
- `[Zvv] = ra.calcZvv(N)`: As in `calcZhh`, but for the vertical polarization.
- `[Zdr] = ra.calcZdr(N)`: Given a drop size distribution `N` (mm^-1 m^-3), estimate the differential reflectivity `Zdr` in units of dB.
- `[rhohv] = ra.calcRhohv(N)`: Given a drop size distribution `N` (mm^-1 m^-3), estimate the co-polar correlation coefficient `rhohv` (unitless, between 0 and 1).
- `[Kdp] = ra.calcKdp(N)`: Given a drop size distribution `N` (mm^-1 m^-3), estimate the specific differential phase `Kdp` in units of deg/km.




