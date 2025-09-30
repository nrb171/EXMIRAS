#! /h/eol/nbarron/miniconda3/envs/pyart-exmiras
# conda install pyart

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
import pyart


folder = "/h/eol/nbarron/work/timrex/spol/20080614/"
file = "cfrad.20080614_110002.000_to_20080614_110716.000_SPOLRVP8_v61_SUR.nc"
radar = pyart.io.read(folder+file)

lon = -16.9
lat = 22.2

col = pyart.util.columnsect.get_field_location(radar, lat, lon)


fig = plt.figure(figsize=(10, 8))




display = pyart.graph.RadarDisplay(radar)
# display.label_yaxis_z()
# display.label_xaxis_time()


display.plot_grid_lines()
# display.plot_ppi(
#     "DBZ",
#     0
# )


display.plot_azimuth_to_rhi(
    "DBZ",
    target_azimuth=314,
    axislabels_flag=True,
    )

# ax.set_xlim(0,35)
# ax.set_ylim(0,2)
display.set_limits(xlim=(0, 60), ylim=(0, 6))
# plt.xlabel("Range (km)", fontsize=12)
# plt.ylabel("Height (km)", fontsize=12)
plt.show()


#! tomorrow:
## one way....
    # figure out how to use these functions in matlab
    
    # grid to high resolution cartesian grid
    # move each level on the grid by a distance based on the wind speed
    # extract simulated RHI
    # use these retrievals to run dsd estimation and EXMIRAS.

## another way...
    # get columns 
    
