import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import my_utils

variable = "Nitrate"  #
depth = 50  # in meter below sealevel
time_frame = "00"  # , "00" for annual. See WOA documentation for other values
resolution = "01"  # see WOA Documentation
data_set = "WOA18"  # see WOA Documentation
color_levels = 11  # number of color levels if z_min_max is not set
z_min_max = [0, 0, 1]  # min, max, step, ignored if min = max
projection = ccrs.Robinson(central_longitude=-150)
fig_size = (10, 7)  # in inches
output_format = "png"
color_scheme = "RdBu_r"
color_scheme = "magma"

# load WOA data
cube, cb_label = my_utils.load_woa_data(
    variable,
    time_frame,
    resolution,
    data_set,
)

# get the nearest depth where data is available
depth_level = my_utils.get_depth_interval(depth)
# extract data at a given level
data = cube[0, depth_level, ...]  # Slice singleton time and first level.

# plot data
fig, ax = my_utils.make_map(
    data,
    projection,
    fig_size,
    cb_label,
    color_scheme,
    z_min_max,
    color_levels,
)

# add title, save figure, and show plot
title = f"{data_set}, {variable}, {depth} mbsl, {time_frame}"
fig_name = f"{data_set}_{variable}_{depth}_{time_frame}_{color_scheme}.{output_format}"
plt.title(title)
fig.savefig(fig_name)
plt.show()
