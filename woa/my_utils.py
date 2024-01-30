def get_depth_interval(depth):
    """Get WOA depth interval from actual depth.  Improve by pre definig step
    size and interval boundaries, rather than setting them manually
    """
    if depth >= 0 and depth <= 125:
        interval = int(round(depth / 5) + 1)
    elif depth <= 500:
        interval = int(round((depth - 125) / 25) + 22)
    elif depth <= 2000:
        interval = int(round((depth - 500) / 50) + 37)
    elif depth <= 5500:
        interval = int(round((depth - 2100) / 100) + 67)
    elif depth > 5000:
        ValueError("WOA is limited to 5500 mbsl")
    else:
        ValueError("This should never happen")

    return interval

def fix_units(cube, field, filename):
    """ WOA data has units as mol_per_kilogram, but it should
    be mol per kilogram
    """
    if "invalid_units" in cube.attributes:
        iu = cube.attributes["invalid_units"]
        cube.units = iu.replace("_", " ")

def load_woa_data(variable, time_frame, resolution, data_set):
    """setup the urls where the data can be found, do some sanity checks on the
    input, download data if necessary, otherwise use local copy, import data,
    and create data cube for a given variable name.

    Parameters:
      variable = name of the data variable, e.g. Oxygen
      depth_level = depth relative to seas surface. See the WOA documentation
      time_frame:str = "00" i.e. annual, monthly etc. see WOA Documentation
      resolution :str = "01" for 1-degree, see WOA Documentation
      data_set :str = "WOA18"
    """
    import urllib.request
    import pathlib as pl
    import iris

    iris.FUTURE.datum_support = True

    res_sep_dict = {"01": "1.00"}
    if resolution in res_sep_dict:
        rsd = res_sep_dict[resolution]
    else:
        ValueError(f"\n no defintion for resolution = {resolution}\n")

    sp_dict = {
        "Phosphate": ["p_an", r"PO$_4$ [$\mu$mol/kg]", "phosphate", "p"],
        "Oxygen": ["o_an", r"O$_2$ [$\mu$mol/kg]", "oxygen", "o"],
        "Nitrate": ["n_an", r"NO$_3$ [$\mu$mol/kg]", "nitrate", "n"],
        "Temperature": ["t_an", r"T [$^{\circ}$ C]", "temperature", "t"],
    }

    if variable in sp_dict:
        var_name = sp_dict[variable][0]
        cb_label = sp_dict[variable][1]
        dir_name = sp_dict[variable][2]
        d_name = sp_dict[variable][3]
        if variable != "Temperature":
            subdir = "all"
        else:
            subdir = "A5B7"  # 2005 - 2017 data
        url = (
            f"https://www.ncei.noaa.gov/data/oceans/woa/{data_set}/"
            f"DATA/{dir_name}/netcdf/{subdir}/{rsd}/"
            f"{data_set.lower()}_{subdir}_{d_name}{time_frame}_{resolution}.nc"
        )
    else:
        raise ValueError(f"\n Unknown variable name {variable}\n")

    print(url)

    fn: str = url.split("/")[-1]  # file name
    cwd: pl.Path = pl.Path.cwd()  # get the current working directory
    fqfn: pl.Path = pl.Path(f"{cwd}/{fn}")  # fully qualified file name

    if not fqfn.exists():  # download data if necessary
        print("Downloading data, this may take a while")
        urllib.request.urlretrieve(url, fqfn)
        print("Downloading finished")
        if not fqfn.exists():  # verify download
            raise FileNotFoundError(f"Cannot find file {fqfn}")

    # import data
    print("Importing Data")
    cubes = iris.load(
        str(fqfn), callback=fix_units
    )  # iris cant handle the path object?

    # find the cube id that matches the desired data
    for k, c in enumerate(cubes):
        if c.var_name == var_name:
            cube_id = k

    cube = cubes[cube_id]  # extract data
    print("Data imported")
    return cube, cb_label

def make_map(cube, projection, figsize, cb_label, color_scheme, z_min_max, levels):
    """Set up map parameters"""
    import iris.plot as iplt
    import matplotlib.cm as mpl_cm
    import matplotlib.pyplot as plt
    import cartopy.feature as cfeature
    import numpy as np

    # create figure instance
    fig, ax = plt.subplots(
        figsize=figsize, subplot_kw=dict(projection=projection), dpi=120
    )
    # color mapping
    color_map = mpl_cm.get_cmap(color_scheme)
    # make plot
    if z_min_max[0] == z_min_max[1]:
        cs = iplt.contourf(cube, cmap=color_map, levels=levels, extend="both")
    else:
        levels = np.arange(z_min_max[0], z_min_max[1] + z_min_max[2], z_min_max[2])
        cs = iplt.contourf(
            cube,
            cmap=color_map,
            vmin=z_min_max[0],
            vmax=z_min_max[1],
            levels=levels,
            extend="both",
        )

    ax.coastlines()  # add coastlines
    ax.add_feature(cfeature.LAND, facecolor="0.75")  # mask land area
    ax.gridlines()

    # add color bar legend
    cb = fig.colorbar(
        cs,
        extend="both",
        shrink=0.5,
        pad=0.02,
        orientation="horizontal",
        fraction=0.1,
    )

    cb.ax.set_xlabel(cb_label)

    return fig, ax
