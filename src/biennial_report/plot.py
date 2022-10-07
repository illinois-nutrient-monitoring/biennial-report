from matplotlib.offsetbox import AnchoredText

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib

def percentage_scale(baseline, ax, diff=False):
    lim = ax.get_ylim()
    ax2 = ax.twinx()
    delta = (lim - baseline)/baseline * 100
    if diff:
        delta = delta + 100
    ax2.set_ylabel('[percentage]')
    ax2.set_ylim(delta)
    return ax2

def annotate_percentage_difference(recent_mean, baseline_mean, loc, ax):
    if loc is None:
        loc='upper right'

    delta = (recent_mean - baseline_mean)/ baseline_mean * 100
    annotate_text = "Change: {:.1f}%".format(delta.values)
    #at = AnchoredText(annotate_text, prop=dict(size=10), frameon=False, loc=loc)
    at = AnchoredText(annotate_text, frameon=False, loc=loc)
    ax.add_artist(at)



def print_period(years):
    return f'{years[0]}-{years[-1]}'

def title_axes(ax):
    ax.set_xlabel( ax.get_xlabel().title())
    ax.set_ylabel( ax.get_ylabel().title())


def annotate_difference(delta, loc, ax):
    if loc is None:
        loc='upper right'

    #delta = (recent_mean - baseline_mean)/ baseline_mean * 100
    annotate_text = "Change: {:.1f}%".format(delta.values)
    #at = AnchoredText(annotate_text, prop=dict(size=10), frameon=False, loc=loc)
    at = AnchoredText(annotate_text, frameon=False, loc=loc)
    ax.add_artist(at)


def mean_and_se(ds, dim):
    mean = ds.mean(dim=dim)
    se = ds.std(dim=dim) / np.sqrt(ds.sizes[dim])
    return mean, se


def plot_mean_error(ds, period, color, ax, linewidth=2, extend=None):
    # ds.sel(river=labels(network2)).sum(dim='river').sel(year=ds.year.dt.year.isin(period))
    mean, se = mean_and_se(ds, dim='year')

    xmin = pd.to_datetime(str(period[0]))
    xmax = pd.to_datetime(str(period[-1]))
    ax.hlines(y=mean, xmin=xmin, xmax=xmax, linewidth=linewidth, color=color, label=print_period(period))
    ax.hlines(y=mean+1.96*se, xmin=xmin, xmax=xmax, linewidth=1, ls='-', color=color)
    ax.hlines(y=mean-1.96*se, xmin=xmin, xmax=xmax, linewidth=1, ls='-', color=color)

    if extend:
        extend = pd.to_datetime(str(extend))
        ax.hlines(y=mean, xmin=xmax, xmax=extend, linewidth=linewidth, ls='dotted', color=color)


def running_average_plot(ds1,
                         period1,
                         period2,
                         ds2=None,
                         ax=None,
                         loc=None):
    if ax is None:
        fig, ax = plt.subplots()

    p1_ds = ds1.sel(year=ds1.year.dt.year.isin(period1))
    p1_mean = p1_ds.mean()

    if ds2 is None:
        temp = ds1
    else:
        temp = ds2

    p2_ds = temp.sel(year=temp.year.dt.year.isin(period2))
    p2_mean = p2_ds.mean()

    
    plot_mean_error(ds=p1_ds,
                    period=period1,
                    color='red',
                    extend=period2[0],
                    ax=ax)


    plot_mean_error(ds=p2_ds,
                    period=period2,
                    color='black',
                    ax=ax)

    ds1.plot(marker='o', color='k',
             linestyle='None',
             markeredgecolor='k',
             markerfacecolor='None',
             ax=ax,
             label=ds1.attrs.get('label'))

    if ds2 is not None:
        ds2.plot(marker='o', color='k',
                 linestyle='None',
                 markeredgecolor='k',
                 markerfacecolor='k',
                 ax=ax,
                 label=ds2.attrs.get('label'))
    
    annotate_percentage_difference(p2_mean, p1_mean, loc, ax)
    title_axes(ax)


def append_total(ds, label='Statewide'):
    total = ds.sum(dim='river').assign_coords({"river":label})
    return xr.concat([total, ds], dim='river')


def plot_by_basin(period_ds,
                  ax=None, color='k',
                  compute_total=True,
                  scale_total=1,
                  mode='load',
                  da=None):
    
    if ax is None:
        fig, ax = plt.subplots()
        
    if compute_total:
        period_ds = append_total(period_ds, scale_total=scale_total, mode=mode, da=da)
    
    title = period_ds.name
    mean, se = mean_and_se(period_ds, dim='year')
    mean = mean.to_series()
    se = se.to_series()
    
    mean.plot.bar(ax=ax, yerr=se*1.96, color=color, edgecolor='k', capsize=4)

    units = period_ds.pint.units
    ax.set_ylabel(f'{title} [{units}]'.capitalize())
    ax.set_xlabel('')
    ax.axhline(y=0, lw=1, color='k')    

def plot_change_by_basin(period1_ds, period2_ds,
                         ax=None, color='k',
                         compute_total=True,
                         scale_total=1,
                         mode='load',
                         da=None):
    '''
    ds1 = ambient_loads[parameter].sel(year=ambient_loads.year.dt.year.isin(baseline_years)).mean(dim='year')
    ds2 = supergage_loads[parameter].sel(year=supergage_loads.year.dt.year.isin(study_years)).mean(dim='year')
    '''
    if ax is None:
        fig, ax = plt.subplots()

    title = period1_ds.name
    

    # append total
    if compute_total:
        period1_ds = append_total(period1_ds, scale_total=scale_total, mode=mode, da=da)
        period2_ds = append_total(period2_ds, scale_total=scale_total, mode=mode, da=da)

    mean1, se1 = mean_and_se(period1_ds, dim='year')
    mean2, se2 = mean_and_se(period2_ds, dim='year')

    difference = (mean2 - mean1).to_series()
    difference_se = np.sqrt(se1**2 + se2**2).to_series()


    #difference.plot.bar(ax=ax, color=color, edgecolor='k', capsize=4)
    difference.plot.bar(ax=ax, yerr=difference_se*1.96, color=color, edgecolor='k', capsize=4)

    units = period1_ds.pint.units
    ax.set_ylabel(f'{title} [{units}]'.capitalize())
    ax.set_xlabel('')
    ax.axhline(y=0, lw=1, color='k')
    #format_axes(ax)

    
def append_total(ds, label='Statewide', scale_total=1, mode='load', da=None):
    #XXXchange to warn
    #assert not (scale_total!=1 and mode=='yield'), "Don't use scale_total with yield=True"
    
    if mode=='load':
        total = ds.sum(dim='river') * scale_total

    elif mode == 'yield':
        total = ds.sum(dim='river')/da.sum()
        ds = ds/da

    return xr.concat([total.assign_coords({"river":label}), ds], dim='river')

def format_yaxis(ax):
    ax.get_yaxis().set_major_formatter(
    matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ','))) 

def format_xaxis(ax):
    ax.get_xaxis().set_major_formatter(
    matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ','))) 

    
#def plot_change_by_basin(ds,
#                         ax=None, color='k',
#                         compute_total=True,
#                         mode='load',
#                         da=None):
#    '''
#    ds1 = ambient_loads[parameter].sel(year=ambient_loads.year.dt.year.isin(baseline_years)).mean(dim='year')
#    ds2 = supergage_loads[parameter].sel(year=supergage_loads.year.dt.year.isin(study_years)).mean(dim='year')
#    '''
#    if ax is None:
#        fig, ax = plt.subplots()
#
#    title = ds.name
#
#    # append total
#    if compute_total:
#        ds = append_total(ds, mode=mode, da=da)
#        #period2_ds = append_total(period2_ds, mode=mode, da=da)
#
#    #mean, se = mean_and_se(ds, dim='year')
#    #mean2, se2 = mean_and_se(period2_ds, dim='year')
#    dim='year'
#    mean = ds.mean(dim=dim)
#    se = ds.std(dim=dim) / np.sqrt(ds.sizes[dim])
#    #difference = (mean2 - mean1).to_series()
#    #difference_se = np.sqrt(se1**2 + se2**2).to_series()
#    mean = mean.to_series()
#    se = se.to_series()
#
#    mean.plot.bar(ax=ax, yerr=se*1.96, color=color, edgecolor='k', capsize=4)
#
#    units = ds.pint.units
#    ax.set_ylabel(f'{title} [{units}]'.capitalize())
#    ax.set_xlabel('')
#    ax.axhline(y=0, lw=1, color='k')
#    #format_axes(ax)
