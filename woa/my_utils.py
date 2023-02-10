def load_woa_data(species, depth_level):
    """setup the urls where the data can be found, do some sanity checks on the
    input, download data if necessary, otherwise use local copy, import data,
    and create data cube for a given variable name.

    Parameters:
      species = name of the data variable, e.g. Oxygen
      depth_level = depth relative to seas surface. See the WOA documentation 
    """
    import urllib.request
    import pathlib as pl
    import iris

    iris.FUTURE.datum_support=True
    
    ou_url = "https://www.ncei.noaa.gov/data/oceans/woa/WOA18/DATA/o2sat/netcdf/all/1.00/woa18_all_O00_01.nc"

    p_url = "https://www.ncei.noaa.gov/data/oceans/woa/WOA18/DATA/phosphate/netcdf/all/1.00/woa18_all_p00_01.nc"

    o_url = "https://www.ncei.noaa.gov/data/oceans/woa/WOA18/DATA/oxygen/netcdf/all/1.00/woa18_all_o00_01.nc"

    T_url = "https://www.ncei.noaa.gov/data/oceans/woa/WOA18/DATA/temperature/netcdf/decav/1.00/woa18_decav_t00_01.nc"

    n_url = "https://www.ncei.noaa.gov/data/oceans/woa/WOA18/DATA/nitrate/netcdf/all/1.00/woa18_all_n00_01.nc"

    sp_dict = {
        "Phosphate": ["p_an", r"PO$_4$ [$\mu$mol/kg]", p_url],
        "Oxygen": ["o_an", r"O$_2$ [$\mu$mol/kg]", o_url],
        "Oxygen_utilisation": ["O_an", r"O$_2$ [$\mu$mol/kg]", o_url],
        "Nitrate": ["n_an", r"NO$_3$ [$\mu$mol/kg]", n_url],
        "Temperature": ["t_an", r"T [$^{\circ}$ C]", T_url],
    }

    if species in sp_dict:
        var_name = sp_dict[species][0]
        cb_label = sp_dict[species][1]
        url = sp_dict[species][2]
    else:
        raise ValueError(f"\n Unknown variable name {species}\n")

    fn: str = url.split("/")[-1]  # file name
    cwd: pl.Path = pl.Path.cwd()  # get the current working directory
    fqfn: pl.Path = pl.Path(f"{cwd}/{fn}")  # fully qualified file name

    if not fqfn.exists():  # download data if necessary
        urllib.request.urlretrieve(url, fqfn)
        if not fqfn.exists():  # verify download
            raise FileNotFoundError(f"Cannot find file {fqfn}")

    # import data
    cubes = iris.load(str(fqfn))  # iris cant handle the pathe object

    # find the cube id that matches the desired data
    for k, c in enumerate(cubes):
        if c.var_name == var_name:
            cube_id = k

    cube = cubes[cube_id]  # extract data
    c = cube[depth_level, 0, ...]  # Slice singleton time and first level.
    print("\n Data imported\n")
    return c, cb_label

def make_map(cube, projection, figsize, cb_label, color_scheme):
    """Set up map parameters"""
    import iris.plot as iplt
    import matplotlib.cm as mpl_cm
    import matplotlib.pyplot as plt
    import cartopy.feature as cfeature

    # create figure instance
    fig, ax = plt.subplots(
        figsize=figsize, subplot_kw=dict(projection=projection), dpi=120
    )
    # plot data
    color_map = mpl_cm.get_cmap(color_scheme)
    cs = iplt.pcolormesh(cube, cmap=color_map)

    ax.coastlines()  # add coastlines
    ax.add_feature(cfeature.LAND, facecolor="0.75")  # mask land area

    # add color bar legend
    cbar = dict(
        extend="both", shrink=0.5, pad=0.02, orientation="horizontal", fraction=0.1
    )
    cb = fig.colorbar(cs, **cbar)
    cb.ax.set_xlabel(cb_label)

    return fig, ax
