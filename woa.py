import os
# os.environ['PROJ_LIB'] = '/path/to/env/share/proj'
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import my_utils

species = "Oxygen"  #
depth_level = 0  # see documentation
time_frame = "annual"  # Match with your download
data_set = "WOA18"  # Match with your download
fig_size = (10, 7)  # in inches
color_scheme = "brewer_RdBu_11"

# load data
data, cb_label = my_utils.load_woa_data(species, depth_level)
# plot data
fig, ax = my_utils.make_map(data, ccrs.Robinson(), fig_size, cb_label, color_scheme)
# add title, save figure, and show plot
title = f"{data_set} {species} {depth_level} {time_frame}"
fig_name = f"{data_set}_{species}_{depth_level}_{time_frame}_{color_scheme}.png"
plt.title(title)
fig.savefig(fig_name)
plt.show()
