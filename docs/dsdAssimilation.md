# `dsdAssimilation.m`
This document describes the `dsdAssimilation` class, which is used to assimilate drop size distribution (DSD) data into a microphysical model based on dual-polarization radar observations. `da` refers to the `dsdAssimilation` object. The `exmiras` class has hooks into this class, but it can be used independently as well.

To initialize a `dsdAssimilation` object, use:
```matlab
da = dsdAssimilation(bandName);
```
## `dsdAssimilation` properties. 
The format for the following properties is `ex.parameterName = defaultValue; %units, size, required`. If there is no `defaultValue`, and there is a `required` flag, then you **must** set the value before the simulation can be run. Bolded properties are required. Italicized properties or `constant` should not be changed unless you know what you are doing. `dependent` properties are calculated from the other properties and should not be set manually.

The properties are grouped into several categories:

### grids and simulation parameters
- `da.zgrid = 25:50:2025 % meters, [sz]`: The vertical grid for the simulation.
- `da.minuteIdx = 7 % integer, [1]`: Use the dsd forwards model LUTs for this minute index. This is the number of minutes after the start of the simulation. Valid values are [4,5,7,10,15].

### observed properties
- `da.RHProfileObs = % seconds, required, [sz]`: Relative humidity profile from observations, interpolated to `zgrid`.
- `da.bandName = 'S' % string`: The radar band for the observations. Options are 'S', 'C', 'X'.
## `dsdAssimilation` methods
- `da.estimateSingleRadarProfile(N0,mu, lambda)`: given a DSD described by the parameters N0, mu, and lambda, estimate the radar reflectivity and differential reflectivity profiles using the forward model LUTs. This takes into account the relative humidity profile.

    Usage:
    - [dNp] = da.estimateSingleRadarProfile(N0,mu, lambda): returns the drop size distribution profile (mm^-1 m^-3)
    - [_,Zhhp] = da.estimateSingleRadarProfile(N0,mu, lambda): returns dNp and the horizontal reflectivity profile (dBZ)
    - [_,_,Zdrp] = da.estimateSingleRadarProfile(N0,mu, lambda): returns dNp, Zhhp, and the differential reflectivity profile (dB)
    - [_,_,_,Kdp] = da.estimateSingleRadarProfile(N0,mu, lambda): returns dNp, Zhhp, Zdrp, and the specific differential phase profile (deg/km)

- `da.profileOptimizer(ZhhProfileObs, ZdrProfileObs)`: Given observed radar reflectivity and differential reflectivity profiles, find the DSD parameters (N0, mu, lambda) at the top of `da.zgrid` that best fit the observations using the dsd forward model LUTs linearization method. `ZhhProfileObs` and `ZdrProfileObs` should be interpolated to `da.zgrid`. This method uses a cost function that minimizes the difference between observed and simulated radar variables.

    Usage: 
    - [N0opt, muopt, lambdaopt] = da.profileOptimizer(ZhhProfileObs, ZdrProfileObs): returns the optimized DSD parameters.
    - [N0opt, muopt, lambdaopt, fv] = da.profileOptimizer(...): also returns the value of the cost function at the optimum.
    - [N0opt, muopt, lambdaopt, fv, dNopt] = da.profileOptimizer(...): also returns the optimized drop size distribution profile (mm^-1 m^-3).
    - [...] = da.profileOptimizer(..., 'KdpProfileObs', KdpProfileObs): also uses observed specific differential phase profile (deg/km) in the cost function. This has mixed success in improving the optimization results.

   Usage:
    - [N0opt, muopt, lambdaopt] = da.profileOptimizer(ZhhProfileObs, ZdrProfileObs, N0guess, muguess, lambdaguess): returns the optimized DSD parameters.
    - [N0opt, muopt, lambdaopt, dNopt] = da.profileOptimizer(...): also returns the optimized drop size distribution profile (mm^-1 m^-3).

- `da.pointOptimizer(ZhhPointObs, ZdrPointObs)`: Given observed radar reflectivity and differential reflectivity at a single height, find the DSD parameters (N0, mu, lambda) that best fit the observations. This does not use the forward model LUTs linearization method, and is closer to the methodology discussed in Zhang et al. (2001).
    Zhang, G., J. Vivekanandan, and E. Brandes, 2001: A method for estimating rain rate and drop size distribution from polarimetric radar measurements. IEEE Trans. Geosci. Remote Sensing, 39, 830–841, https://doi.org/10.1109/36.917906.   

    Usage:
    - [N] = da.pointOptimizer(ZhhPointObs, ZdrPointObs): returns the optimized DSD (a 1 x numel(D) array).
    - [N, N0opt, muopt, lambdaopt] = da.pointOptimizer(...): also returns the optimized DSD parameters. 

## Conversion between radar variables/DSD variables to DSD parameters
- `[N,N0,mu,gamma]=da.getNFromZhhZdr(Zhh, Zdr)`: Given radar reflectivity (dBZ) and differential reflectivity (dB), return the DSD parameters (N0, mu, lambda) using the method of Zhang et al. (2001).
    - note: `da.pointOptimizer` is just an alias of this function.
- `[N,N0, mu, gamma] = getNFromZhhZdrDm(Zhh, Zdr, Dm)`: Given radar reflectivity (dBZ), differential reflectivity (dB), and mean diameter (mm), return the DSD parameters (N0, mu, lambda).
- `N0 = getN0FromZhhMuLambda(obj, Zhh, mu, lambda)`: Given radar reflectivity (dBZ), mu, and lambda, return N0. This method is helpful for converting a DSD of one reflectivity to another without changing mu and lambda.
- `[N0,mu,lambda] = getN0MuLambdaFromN(obj, N)`: Given a DSD (mm^-1 m^-3), return the DSD parameters (N0, mu, lambda).

## Using the `dsdAssimilation` class
To use the `dsdAssimilation` class, you typically follow these steps:
1. Initialize the `dsdAssimilation` object with the desired radar band.
2. Set the required properties, such as `RHProfileObs`.
3. Use the `profileOptimizer` or `pointOptimizer` methods to estimate the DSD parameters from observed radar data.


