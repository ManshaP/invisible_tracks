
from pyhdf.SD import SD, SDC


def delete_no_overlap(filelist, extent):
    delfiles=[]
    for filename in filelist:
        #print(filename)
        try: dataset = SD(filename, SDC.READ)
        except: 
            delfiles.append(filename)
            continue
        lats = dataset.select('Latitude').get()
        lons = dataset.select('Longitude').get()
        #print(lats.min(), lats.max(), lons.min(), lons.max())
        if (lats.max() < extent[0][0]  or 
            lats.min() > extent[0][1] or 
            lons.max() < extent[1][0] or 
            lons.min() > extent[1][1] or 
            lons.max()>170 ):
            #print(filename)
            delfiles.append(filename)
        #else: print(filename)
    
    for filename in delfiles:
        filelist.remove(filename)
    return filelist


def UngriddedData_from_data_frame(df, cols, names=None, air_pressure_units=None):
    """
    Create an UngriddedData object from the cols of a datafame (df)

    :param df: The input dataframe
    :param cols: The columns to extract (note that coordinates are dealt with automatically)
    :param names: The names to give the date objects (if different to the column names)
    :param air_pressure_units: Optional air pressure units which the data was ORIGINALLY in. The output is always hPa
    :return UngriddedDataList: List of UngriddedData objects, one for each column specified
    """
    from cis.data_io.ungridded_data import UngriddedData, UngriddedDataList, Metadata
    from cis.data_io.Coord import Coord, CoordList
    from cis.data_io.write_netcdf import types as valid_types
    from cis.time_util import cis_standard_time_unit
    from cf_units import Unit
    from numpy import ma

    fill_value = -9e30
    # define our function to perform the case-insensitive search
    def find_col_name(name):
        col_list = list(df)
        try:
            # this uses a generator to find the index if it matches, will raise an exception if not found
            return col_list[next(i for i, v in enumerate(col_list) if v.lower() == name)]
        except:
            return ''

    lat_col_name = find_col_name('latitude')
    lon_col_name = find_col_name('longitude')
    alt_col_name = find_col_name('altitude')
    pres_col_name = find_col_name('air_pressure')

    coords = CoordList()
    out_data = UngriddedDataList()
    numpy_time_unit = Unit('ns since 1970-01-01T00:00Z')
    time_vals = numpy_time_unit.convert(df.index.values.astype('float64'), cis_standard_time_unit)
    coords.append(Coord(time_vals, Metadata(standard_name='time', units=str(cis_standard_time_unit))))
    coords.append(Coord(df[lat_col_name].values, Metadata(standard_name='latitude', units='degrees_north')))
    coords.append(Coord(df[lon_col_name].values, Metadata(standard_name='longitude', units='degrees_east')))
    df = df.drop([lat_col_name, lon_col_name], axis=1)
    if alt_col_name in df:
        coords.append(Coord(df[alt_col_name].values, Metadata(standard_name='altitude', units='meters')))
        df = df.drop([alt_col_name], axis=1)
    if pres_col_name in df:
        air_pressure_units = air_pressure_units if air_pressure_units is not None else 'hPa'
        pres_data = Unit(air_pressure_units).convert(df[pres_col_name].values, 'hPa')
        coords.append(Coord(pres_data, Metadata(standard_name='air_pressure', units='hPa')))
        df = df.drop([pres_col_name], axis=1)

    # Check the names and cols match up if present
    if (cols and names) and (len(cols) != len(names)):
        raise ValueError()
    elif not names:
        names = cols

    for col, _name in zip(cols, names):
        if str(df[col].values.dtype) in valid_types.keys():
            if df[col].isnull().any():
                # Note we specify the mask explitly for each column because we could have coiord values which are valid
                # for some data variables and not others
                data = ma.array(df[col].values, mask=df[col].isnull(), fill_value=fill_value)
                data[data.mask] = fill_value
            else:
                data = df[col].values
            meta = Metadata(long_name=col, name=_name)
            if str(df[col].values.dtype) != 'object':
                meta.missing_value = fill_value
            out_data.append(UngriddedData(data, meta, coords.copy()))
    return out_data
