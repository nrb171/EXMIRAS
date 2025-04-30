Table of Contents:
- [EXMIRAS](#exmiras)
  - [Installation](#installation)
    - [Dependencies](#dependencies)
- [Instructions](#instructions)
- [Running the test cases](../docs/runningTheTestCases.md)
  
# EXMIRAS
EXMIRAS is a MATLAB-based software package for simulating liquid rainfall and dual-pol radar observations of that rainfall. 

## Installation
Clone the repository on a machine with MATLAB installed: 
```
git clone https://github.com/nrb171/EXMIRAS
```

### Dependencies
- The main scripts, [`../scripts/exmiras.m`](../scripts/exmiras.m) and [`../scripts/thermo.m`](../scripts/thermo.m) only require MATLAB. These are the core of the EXMIRAS package. 
- There are several optional python scripts with various dependencies. Please refer to the installation instructions for each package below.
    - I recommend building a virtual environment for each of the scripts below, since [`pytmatrix`](https://github.com/jleinonen/pytmatrix) is an older package that may conflict with other packages.
    - [`./scripts/helpers/tmatrixInterface.py`](../scripts/helpers/tmatrixInterface.py) requires the [`pytmatrix`](https://github.com/jleinonen/pytmatrix), and [`numpy`](https://pypi.org/project/numpy/). This script was used to compute the scattering properties of the hydrometeors in the simulations and generates the LUTs in [`../data/LUTs`](../data/LUTs/). You should not need to regenerate these LUTs, but the script is included for reference.
    - [`../scripts/helpers/processASOS.py`](../scripts/helpers/processASOS.py) requires [`pandas`](https://pandas.pydata.org/), [`numpy`](https://pypi.org/project/numpy/), and [`metpy`](https://pypi.org/project/MetPy/). This script was used to process the metar data for the case study in the paper. If you are using a different dataset, you will need to modify the script accordingly.

## Instructions
From this point forward, `exmiras` refers to the `exmiras.m` class, while `EXMIRAS` refers to the entire package.
1. First, create a new script and call the `exmiras` class:
    ```matlab
    ex = exmiras;
    ```
### `exmiras` properties. 
2. Before running the simulation, you will need to set many of the properties for the simulation. You can view the list of all properties by looking at the top of the [`exmiras.m`](../scripts/exmiras.m) or [`thermo.m`](../scripts/thermo.m) scripts. 

    The format for the following properties is `ex.parameterName = defaultValue; %units, size, required`. If there is no `defaultValue`, and there is a `required` flag, then you **must** set the value before the simulation can be run. Bolded properties are required. Italicized properties or `constant` should not be changed unless you know what you are doing. `dependent` properties are calculated from the other properties and should not be set manually.

    The properties are grouped into several categories:

    #### simulation properties
    - `ex.dt = 10 % seconds`: The time step for the simulation. All non-fixed properties are updated at this time step.
    - `ex.st % seconds`: The simulation time. This is the total time for the simulation.
    - **`ex.xgrid % meters, required, [xsize]`**: The x grid for the simulation. 
    - **`ex.ygrid % meters, required, [ysize]`**: The y grid for the simulation.
    - **`ex.zgrid % meters, required, [zsize]`**: The z (vertical) grid for the simulation.
    - *`ex.nBins = 250 % integer, constant`*: The number of droplet size bins for the simulation.  
    Note: the domain will be set as `[sx, sy, sz, nbins]` where sx/sy/sz are the numel(sizex/y/z).

    #### number concentration properties
    These properties control the drop size distribution (DSD) for the simulation. It is recommended to set these values automatically using the `ex = ex.initFromReflectivity(dBZi, Zdri);` class method where `dBZi` is the reflectivity in dBZ and `Zdri` is the differential reflectivity in dBZ. However, these can all be set manually if desired.
    - **`ex.N % m^-3 mm^-1, [sx, sy, sz, nbins], required`**: The number concentration for the simulation. This is the number of droplets per cubic meter per millimeter of diameter.
    - **`ex.N0 % m^-3 mm^-1, [sx, sy, sz, nbins], required`**: the dsd from the previous time step. This is used to calculate several rates of change in the simulation.
    - *`ex.N00 = ex.N; % m^-3 mm^-1, [sx, sy, sz, nbins], required`*: The initial dsd for the simulation. This is used to calculate several rates of change in the simulation.
    - **`ex.mu % shape parameter, required`**: The shape parameter for the gamma distribution. This is used to calculate the drop size distribution (DSD) for the simulation.
    - **`ex.gamma % scale parameter, required`**: The scale parameter for the gamma distribution. 
    - **`ex.NP % m^-3 mm^-1, [sx, sy, sz, nbins], required`**: This is the number concentration from the 'parent' field. This sets the initial number concentration for the simulation and is replenished at the end of each time step. This essentially is the source of the rainfall. <u>**This value is not updated from the `ex.initFromReflectivity` method.**</u>

    #### radar properties
    These properties store the radar simulations or the radar settings.
    - `ex.lambda = 111 % mm`: The wavelength for the radar simulation. This can also be set using the name of the wavelength from the `ex = ex.initFromLambdaName(lambdaName);` class method. The options/wavelengths are as follows:
        - `bandName = ["S",    "C",    "X",        "Ku",       "Ka"       ];`
        - `wavelength = [111,  53.5,   33.3,       22,         8.43       ];`
    - *`ex.dpp % constant`*: the structure holding the scattering properties from the LUTs. 
    - *`ex.Zhh % dBZ, [sx, sy, sz], dependent`*: the horizontally-polarized reflectivity.
    - *`ex.Zvv % dBZ, [sx, sy, sz], dependent`*: the vertically-polarized reflectivity.
    - *`ex.Zdr % dBZ, [sx, sy, sz], dependent`*: the differential reflectivity.
    - *`ex.kdp % deg/km, [sx, sy, sz], dependent`*: the specific differential phase.
    - *`ex.rhohv %correlation, [sx, sy, sz], dependent`*: the cross-correlation coefficient.
    - *`ex.RR % mm/hr, [sx, sy, sz], dependent`*: the rain rate.
    
    #### droplet properties
    - *`ex.D % mm, [nbins], dependent`*: the characteristic droplet diameter in each bin.
    - *`ex.Dw % mm, [nbins], dependent`*: the width of each droplet size bin.
    - *`ex.Dmax % mm, [1], constant`*: the maximum droplet diameter.
    - *`ex.M %kg , [nbins], dependent`*: the mass of each droplet in each bin.
    - *`ex.vt % m/s, [nbins], dependent`*: the terminal velocity of each droplet in each bin.

    #### atmospheric thermodynamic and kinematic properties
    - **`ex.T % K, [sx, sy, sz], required`**: the initial temperature for the simulation. Definition in `thermo.m`
    - **`ex.p % hPa, [sx, sy, sz], required`**: the initial pressure for the simulation. Definition in `thermo.m`
    - **`ex.pv % hPa, [sx, sy, sz], required`**: the initial partial pressure of water vapor.
    - *`ex.rhoa % kg/m^3, [sx, sy, sz], dependent`*: the density of air. Definition in `thermo.m`
    - *`ex.m % kg/m^3, [sx, sy, sz], dependent`*: total mass of liquid water in the air volume.
    - *`ex.mv % kg/m^3, [sx, sy, sz], dependent`*: total mass of water vapor in the air volume.
    - *`ex.theta % K, [sx, sy, sz], dependent`*: the potential temperature.
    - **`ex.u/v/w % m/s, [sx, sy, sz], required`**: the horizontal and vertical wind speed. Only the vertical wind speed is used in the simulation, u/v may be used in the future.
    
    #### rates of change
    - *`ex.dmevap % kg/m^3, [sx, sy, sz, nbins], dependent`*: the change in mass due to evaporation.
    - *`ex.dDevap % mm, [sx, sy, sz, nbins], dependent`*: the change in droplet size due to evaporation.
    - *`ex.dNevap % m^-3 mm^-1, [sx, sy, sz, nbins], dependent`*: the change in number concentration due to evaporation.
    - *`ex.dNadvect % m^-3 mm^-1, [sx, sy, sz, nbins], dependent`*: the change in number concentration due to advection. <u>**this method is not yet implemented.**</u>
    - *`ex.dNfallout % m^-3 mm^-1, [sx, sy, sz, nbins], dependent`*: the change in number concentration due to fallout.
    - *`ex.dNcoal % m^-3 mm^-1, [sx, sy, sz, nbins], dependent`*: the change in number concentration due to collision-coalescence.
    - *`ex.dTevap % K, [sx, sy, sz], dependent`*: the change in atmospheric temperature from evaporative processes.

### running the simulation
3. Once all of the properties have been set, you can run the simulation using the `ex = ex.integrate();` class method. This will simulate a single time step and update all of the properties at each time step. If you want to record the simulation's result at each time steps, using something like the following code snippet:

    ```matlab
    T = zeros(numSteps, sz);
    p = zeros(numSteps, sz);
    qv = zeros(numSteps, sz);
    pv = zeros(numSteps, sz);
    dBZhh = zeros(numSteps, sz);
    dBZvv = zeros(numSteps, sz);
    theta = zeros(numSteps, sz);
    RR = zeros(numSteps, sz);
    rhohv = zeros(numSteps, sz);
    kdp = zeros(numSteps, sz);
    for i = 1:numSteps
        ev = ex.integrate;
        T(i,:) = squeeze(ex.T(1,1,:))-squeeze(temp(zm(1,1,:)));
        theta(i,:) = ex.theta(1,1,:);
        p(i,:) = ex.p(1,1,:);
        qv(i,:) = ex.qv(1,1,:);
        pv(i,:) = ex.pv(1,1,:);
        dBZhh(i,:) = ex.Zhh(1,1,:);
        dBZvv(i,:) = ex.Zvv(1,1,:);
        RR(i,:) = ex.RR(1,1,:);
        rhohv(i,:) = ex.rhohv(1,1,:);
        kdp(i,:) = ex.kdp(1,1,:);
    end
    ```
    This will record height-time profiles of the simulation's results. You can then plot them using [`../scripts/helpers/plotTimeProfile.m`](../scripts/helpers/plotTimeProfile.m).





    




