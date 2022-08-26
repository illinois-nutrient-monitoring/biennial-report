import xarray as xr

def location_label(site):
    return site['river'] + ' ' + site['location']


def compute_river_load(ds, site, subtract_nested=True):

    site_load = ds.sel(site=site['gage_id'])
    site_load = site_load.drop_vars('site')

    upstream_gage = site.get('upstream_gage')
    
    # if nesting enabled, subtract upstream load
    if subtract_nested and upstream_gage:
        upstream_load = ds.sel(site=upstream_gage).drop_vars('site')
        site_load = site_load - upstream_load

    site_load = site_load.assign_coords(coords={'river':location_label(site)})

    site_load = site_load * site['scale_factor']
    return site_load


def compute_network_loads(ds, network, subtract_nested=True):
    rivers = []

    for site in network:
        if site['gage_id'] in ds.site:
            nested = site.get('nested')
            if not nested or subtract_nested is False:
                rivers.append(compute_river_load(ds, site, subtract_nested))

    return xr.concat(rivers, dim='river')


def kg_to_lbs(x):
    return x * 2.20462

def kg_to_kt(x):
    return x * 1.0e-6

def labels(network, include_nested=False):
    labels = []
    for site in network:
        is_nested = site.get('nested')
        
        if is_nested is None:
            labels.append(location_label(site))
        elif include_nested and is_nested is True:
            labels.append(location_label(site))

#        if is_nested is None or is_nested==nested:
#            labels.append(location_label(site))

    return labels

def gages(network, include_nested=False):
    labels = []
    for site in network:
        is_nested = site.get('nested')

        if is_nested is None:
            labels.append(site['gage_id'])
        elif include_nested and is_nested is True:
            labels.append(site['gage_id'])

    return labels