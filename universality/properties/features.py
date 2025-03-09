"""a module housing logic to identify (highly engineered) features from EoS
"""
__author__ = "Reed Essick (reed.essick@gmail.com)"

#-------------------------------------------------

import os

import numpy as np

from argparse import ArgumentParser

### non-standard libraries
from universality.utils import (io, utils)
from universality.plot.utils import plt

### some additional imports 

from scipy.ndimage import gaussian_filter1d
import matplotlib.ticker as ticker 


#-------------------------------------------------

DEFAULT_EOS_COLUMNS = ['baryon_density', 'pressurec2', 'energy_densityc2', 'cs2c2']
DEFAULT_MACRO_COLUMNS = ['M', 'R', 'I', 'Lambda']

#-------------------------------------------------
# basic utilities for finding features before/after a particular point
#-------------------------------------------------

def find_preceeding(x_ref, x, y):
    '''assumes x is ordered from decreasing to increasing'''
    found = x <= x_ref
    if np.sum(found):
        return y[found][-1]
    else:
        raise RuntimeError('no preceeding value found!')

def find_following(x_ref, x, y):
    '''assumes x is ordered from decreasing to increasing'''
    return find_preceeding(-x_ref, -x[::-1], y[::-1])

#------------------------

def find_running_maxima(x):
    """find the indecies of the largest thing that's been seen up until the current index"""
    return find_running_minima(-np.array(x))

def find_running_minima(x):
    """find the indecies of the smallest thing that's been seen up until the current index"""
    mins = []
    best = +np.infty
    for i, X in enumerate(x):
        if X < best:
            mins.append(i)
            best = X
    return mins

#------------------------

def find_maxima(x):
    """find indecies of strict local maxima"""
    return find_minima(-np.array(x))

def find_minima(x):
    """find indecies of strict local minima"""
    mins = []
    old = +np.infty
    for i in range(len(x)-1):
        new = x[i]
        if (old > new) and (new < x[i+1]):
            mins.append(i)
        old = new
    if x[i+1] > x[i]:
        mins.append(i+1)
    return mins

#------------------------

def find_inclusive_maxima(x):
    """find indecies of local maxima, including cases when the value does not change from one sample to the next
if the derivative vanishes (neighboring samples are identical), then only the entry in the sequence of identical entries is returned
    """
    return find_inclusive_minima(-np.array(x))

def find_inclusive_minima(x):
    """find indecies of local minima, including cases when the value does not change from one sample to the next
if the derivative vanishes (neighboring samples are identical), then only the entry in the sequence of identical entries is returned
    """
    mins = []
    old = +np.infty
    for i in range(len(x)-1):
        new = x[i]
        ### this will keep ONLY the earliest occurance in a sequence where the derivative vanishes
        if (old > new) and (new <= x[i+1]):
            mins.append(i)
        old = new
    if x[i+1] > x[i]:
        mins.append(i+1)
    return mins

#-------------------------------------------------
# Features extracted from moment of inertia relations
#-------------------------------------------------

# templates for column names relevant to MoI features
MAX_CS2C2_TEMPLATE = 'max_cs2c2_%s'
RMAX_CS2C2_TEMPLATE = 'running_max_cs2c2_%s'
MIN_CS2C2_TEMPLATE = 'min_cs2c2_%s'

LOW_CS2C2_TEMPLATE = 'smallest_cs2c2_%s'
HGH_CS2C2_TEMPLATE = 'largest_cs2c2_%s'

MAX_ARCTAN_DLNI_DLNM_TEMPLATE = 'max_arctan_dlnI_dlnM_%s'
MIN_ARCTAN_DLNI_DLNM_TEMPLATE = 'min_arctan_dlnI_dlnM_%s'

DEFAULT_FLATTEN_THR = 0.0
DEFAULT_SMOOTHING_WIDTH = 0.01 ### applied in log(rhoc). This seems to be small enough to not affect much
                               ### but (anectodally) large enough to smooth out typical numeric variation
DEFAULT_DIFF_THR = 0.0
DEFAULT_CS2C2_DROP_RATIO = 0.0

#------------------------

### utility functions for identifying MoI features (including basic feature transformation)

def arctan_transform(rhoc, M, I, flatten_thr=DEFAULT_FLATTEN_THR, smoothing_width=DEFAULT_SMOOTHING_WIDTH):
    """compute the special arctan function we use to feed into our feature identification function
    """
    dlnI_drhoc = utils.num_dfdx(rhoc, I) / I
    dlnM_drhoc = utils.num_dfdx(rhoc, M) / M

    if smoothing_width:
        dlnI_drhoc = smooth(np.log(rhoc), dlnI_drhoc, smoothing_width)
        dlnM_drhoc = smooth(np.log(rhoc), dlnM_drhoc, smoothing_width)

    ### regularlize the derivatives so we don't get oscillatory behavior due to numerical precision when there is no change with respect to rhoc
    ### only filter places where we've touched both dM and dI to account for numerical error
    spurious = (np.abs(dlnI_drhoc*rhoc)<flatten_thr) * (np.abs(dlnM_drhoc*rhoc)<flatten_thr)

    dlnI_drhoc[spurious] = 0
    dlnM_drhoc[spurious] = 0

    ans = np.arctan2(dlnI_drhoc, dlnM_drhoc)
    ans[ans>np.pi/2] -= 2*np.pi ### wrap this around so that all the region where dlnM_drhoc < 0 is "contiguous"

    return ans, (spurious, dlnM_drhoc, dlnI_drhoc)

def smooth(x, y, width):
    v = (width)**2
    ans = np.empty_like(y)
    inds = np.arange(len(x))
    for i, S in enumerate(x):
        weights = np.exp(-0.5*(S-x)**2/v) ### the basic Gaussian kernel
        weights /= np.sum(weights)
        ans[i] = np.sum(weights*y)
    return ans

#------------------------

### functions to extract parameters

DEFAULT_MOI_FEATURE_NAME = 'moi_feature'

def data2moi_features(
        rhoc,
        M,
        I,
        baryon_density,
        cs2c2,
        macro_data,
        macro_cols,
        eos_data,
        eos_cols,
        flatten_thr=DEFAULT_FLATTEN_THR,
        smoothing_width=DEFAULT_SMOOTHING_WIDTH,
        diff_thr=DEFAULT_DIFF_THR,
        cs2c2_drop_ratio=DEFAULT_CS2C2_DROP_RATIO,
        verbose=False,
        debug_figname=None, ### this is a path into which we write the debug figure
    ):
    """
    looks for specific behavior in arctan(dlnI/dlnM) that seems to correspond to phase transitions.
    returns the corresponding values from data along with names based on cols
    """
    names = [DEFAULT_MOI_FEATURE_NAME]
    for tmp in [
            MIN_CS2C2_TEMPLATE,
            MAX_CS2C2_TEMPLATE,
            RMAX_CS2C2_TEMPLATE,
            LOW_CS2C2_TEMPLATE,
            HGH_CS2C2_TEMPLATE,
            MAX_ARCTAN_DLNI_DLNM_TEMPLATE,
            MIN_ARCTAN_DLNI_DLNM_TEMPLATE,
        ]:
        names.append(tmp % 'arctan_dlnI_dlnM')
        names += [tmp%col for col in macro_cols]
        names += [tmp%col for col in eos_cols]
    params = []

    if debug_figname: ### set up lists for plotting
        points = []  # individual points
        regions = [] # spans of density

    #---

    ### compute the absolute value of the curvature, which we use as an indicator variable
    arctan_dlnI_dlnM, (spurious, dlnM_drhoc, dlnI_drhoc) = arctan_transform(
        rhoc,
        M,
        I,
        flatten_thr=flatten_thr,
        smoothing_width=smoothing_width,
    )

    ### find the possible end points as local minima of arctan_dlnI_dlnM
    ends = list(find_inclusive_minima(arctan_dlnI_dlnM)[::-1]) ### reverse so the ones with largest rhoc are first

    if debug_figname:
        points += [(rhoc[end], dict(color='k', marker='o')) for end in ends]

    # discard any local minima that are before the first stable branch
    while len(ends):
        end = ends[-1]
        if np.any(dlnM_drhoc[:end] > 0): ### something is stable before this
             break
        ends = ends[:-1] ### truncate this guy

        if debug_figname: # add to list for plotting
            points.append((rhoc[end], dict(color='b', marker='o')))

    # discard any local minima that are in the final unstable branch
    while len(ends):
        end = ends.pop(0)
        if np.any(dlnM_drhoc[end:] > 0): ### not part of the 'final unstable branch'
            ends.insert(0, end)
            break

        if debug_figname: # add to list for plotting
            points.append((rhoc[end], dict(color='r', marker='s')))

    #---

    if ends: ### we have something to do

        ### if the lowest rhoc available is a minimum, we discard it since we can't tell if it is a local minimum or not
        if ends[-1] == 0:
            ends = ends[:-1]

            if debug_figname:
                points.append((rhoc[0], dict(color='g', marker='*')))

        if ends[0] == len(rhoc)-1: ### same thing with the end point
            ends = ends[1:]

            if debug_figname:
                points.append((rhoc[-1], dict(color='g', marker='*')))

        #---

        ### local minima in sound speed
        min_cs2c2 = np.array(find_inclusive_minima(cs2c2))
        min_cs2c2_baryon_density = baryon_density[min_cs2c2]

        ### local maxima in sound speed
        max_cs2c2 = np.array(find_inclusive_maxima(cs2c2))
        max_cs2c2_baryon_density = baryon_density[max_cs2c2]

        ### global maxima in sound speed up to the current point
        rmax_cs2c2 = find_running_maxima(cs2c2)

        if not rmax_cs2c2: ### no qualifying maxima
            if debug_figure:
                for end in ends:
                    points.append((rhoc[end], dict(color='y', marker='v')))

            ends = [] ### this will make us skip all the ends because we can't match them to starting points

        else:
            rmax_cs2c2 = np.array(rmax_cs2c2)
            rmax_cs2c2_baryon_density = baryon_density[rmax_cs2c2]

        #---

        ### iterate through and grab the following associated with each "end"
        Neos = len(eos_cols)
        Nmac = len(macro_cols)

        for end in ends:
            r = rhoc[end]

            if verbose:
                print('        processing end=%d/%d at rhoc=%.6e'%(end, len(rhoc), r))

            #---

            ### min sound speed preceeding "end"
            try:
                ind_min_cs2c2 = find_preceeding(r, min_cs2c2_baryon_density, min_cs2c2)
            except RuntimeError:
                if verbose:
                    print('            WARNING! could not find local minimum in cs2c2 preceeding local minimum in arctan(dlogI/dlogM)')

                if debug_figname:
                    points.append((r, dict(color='orange', marker='h')))

                continue

            min_r = baryon_density[ind_min_cs2c2]
            min_cs2c2_arctan = np.interp(min_r, rhoc, arctan_dlnI_dlnM)

            #---

            ### max sound speed preceeding minimum sound speed
            try:
                ind_max_cs2c2 = find_preceeding(min_r, max_cs2c2_baryon_density, max_cs2c2)
            except RuntimeError:
                if verbose:
                    print('            WARNING! could not find local maximum in cs2c2 preceeding local minimum in cs2c2')

                if debug_figname:
                    points.append((r, dict(color='orange', marker='h')))

                continue

            max_r = baryon_density[ind_max_cs2c2]

            #---

            ### running max sound speed preceeding maximum sound speed

            try:
                ind_rmax_cs2c2 = find_preceeding(max_r, rmax_cs2c2_baryon_density, rmax_cs2c2)
            except RuntimeError:
                if verbose:
                    print('            WARNING! could not find running maximum in cs2c2 preceeding local minimum in cs2c2')

                if debug_figname:
                    points.append((r, dict(color='c', marker='d')))

                continue

            rmax_r = baryon_density[ind_rmax_cs2c2]
            rmax_cs2c2_arctan = np.interp(rmax_r, rhoc, arctan_dlnI_dlnM)

            ### now find the maximum arctan(...) between rmax_r and r
            selected = (rmax_r<=rhoc)*(rhoc<=r)
            max_arctan_r = rhoc[selected][np.argmax(arctan_dlnI_dlnM[selected])]
            ind_max_arctan = np.arange(len(rhoc))[rhoc==max_arctan_r][0]

            selected = (rmax_r<=baryon_density) * (baryon_density<=r)
            low_cs2c2_r = baryon_density[selected][np.argmin(cs2c2[selected])]
            ind_low_cs2c2 = np.arange(len(cs2c2))[low_cs2c2_r==baryon_density][0]

            hgh_cs2c2_r = baryon_density[selected][np.argmax(cs2c2[selected])]
            ind_hgh_cs2c2 = np.arange(len(cs2c2))[hgh_cs2c2_r==baryon_density][0]

            #---

            ### now perform sanity checks to make sure this candidate passes
            ### iterate until we find an acceptable running max_cs2c2

            # select an earlier rmax_cs2c2 if Delta(arctan) is too small
            diff_arctan = arctan_dlnI_dlnM[ind_max_arctan] - arctan_dlnI_dlnM[end]

            # sound speed must drop by a certain fraction
            drop_ratio = cs2c2[ind_rmax_cs2c2] / cs2c2[ind_low_cs2c2]

            while (diff_arctan < diff_thr) or (drop_ratio < cs2c2_drop_ratio):

                # first, find preceeding local max_cs2c2
                try:
                    truth = max_cs2c2_baryon_density < rmax_r ### move to a local max that is before the current running max
                    jnd_max_cs2c2 = find_preceeding(rmax_r, max_cs2c2_baryon_density[truth], max_cs2c2[truth])
                except RuntimeError:
                    if verbose:
                        print('            WARNING! could not find local maximum in cs2c2 preceeding local minimum in cs2c2')
                    if debug_figname:
                        points.append((r, dict(color='orange', marker='h')))

                    break

                # then update rmax_cs2c2 to be before that
                try:
                    ind_rmax_cs2c2 = find_preceeding(baryon_density[jnd_max_cs2c2], rmax_cs2c2_baryon_density, rmax_cs2c2)
                except RuntimeError:
                    if verbose:
                        print('            WARNING! could not find running maximum in cs2c2 preceeding local minimum in cs2c2')

                    if debug_figname:
                        points.append((r, dict(color='c', marker='d')))

                    break

                rmax_r = baryon_density[ind_rmax_cs2c2]
                rmax_cs2c2_arctan = np.interp(rmax_r, rhoc, arctan_dlnI_dlnM)

                ### now find the maximum arctan(...) between rmax_r and r
                selected = (rmax_r<=rhoc) * (rhoc<=r)
                max_arctan_r = rhoc[selected][np.argmax(arctan_dlnI_dlnM[selected])]
                ind_max_arctan = np.arange(len(rhoc))[rhoc==max_arctan_r][0]

                selected = (rmax_r<=baryon_density) * (baryon_density<=r)
                low_cs2c2_r = baryon_density[selected][np.argmin(cs2c2[selected])]
                ind_low_cs2c2 = np.arange(len(cs2c2))[low_cs2c2_r==baryon_density][0]

                hgh_cs2c2_r = baryon_density[selected][np.argmin(cs2c2[selected])]
                ind_hgh_cs2c2 = np.arange(len(cs2c2))[hgh_cs2c2_r==baryon_density][0]

                # update conditionals
                diff_arctan = arctan_dlnI_dlnM[ind_max_arctan] - arctan_dlnI_dlnM[end]
                drop_ratio = cs2c2[ind_rmax_cs2c2] / cs2c2[ind_low_cs2c2]

            if (diff_arctan < diff_thr) or (drop_ratio < cs2c2_drop_ratio): # we exited from the break statement
                continue

            #--- add parameters at identified points (same ordering as when we construct "names")

            datum = []

            # add local min in cs2c2
            datum.append(np.interp(min_r, rhoc, arctan_dlnI_dlnM))
            datum += [np.interp(min_r, rhoc, macro_data[:,i]) for i in range(Nmac)]
            datum += list(eos_data[ind_min_cs2c2])

            # add local max in cs2c2
            datum.append(np.interp(max_r, rhoc, arctan_dlnI_dlnM))
            datum += [np.interp(max_r, rhoc, macro_data[:,i]) for i in range(Nmac)]
            datum += list(eos_data[ind_max_cs2c2])

            # add running max in cs2c2
            datum.append(np.interp(rmax_r, rhoc, arctan_dlnI_dlnM))
            datum += [np.interp(rmax_r, rhoc, macro_data[:,i]) for i in range(Nmac)]
            datum += list(eos_data[ind_rmax_cs2c2])

            # add lowest cs2c2 within feature
            datum.append(np.interp(low_cs2c2_r, rhoc, arctan_dlnI_dlnM))
            datum += [np.interp(low_cs2c2_r, rhoc, macro_data[:,i]) for i in range(Nmac)]
            datum += list(eos_data[ind_low_cs2c2])

            # add highest cs2c2 within feature
            datum.append(np.interp(hgh_cs2c2_r, rhoc, arctan_dlnI_dlnM))
            datum += [np.interp(hgh_cs2c2_r, rhoc, macro_data[:,i]) for i in range(Nmac)]
            datum += list(eos_data[ind_hgh_cs2c2])

            # add parameters at max arctan
            datum.append(arctan_dlnI_dlnM[ind_max_arctan])
            datum += list(macro_data[ind_max_arctan])
            datum += [np.interp(max_arctan_r, baryon_density, eos_data[:,i]) for i in range(Neos)]

            # add parameters at "end" (min arctan)
            datum.append(arctan_dlnI_dlnM[end])
            datum += list(macro_data[end])
            datum += [np.interp(r, baryon_density, eos_data[:,i]) for i in range(Neos)]

            if debug_figname: # plot surviving candidates
#                regions.append(((rmax_r, max_r, min_r, r), dict(color='k', marker='.', linewidth=2, alpha=0.25)))
                regions.append(((rmax_r, max_r, min_r, max_arctan_r, r), dict(marker='.', linewidth=2)))

            #---

            params.append(datum)

    #---

    if debug_figname: ### make a plot

        if verbose:
            print('plotting')

        fig = data2moi_features_figure(
            points,
            regions,
            rhoc,
            M,
            I,
            baryon_density,
            cs2c2,
            arctan_dlnI_dlnM,
            dlnM_drhoc,
            dlnI_drhoc,
            '%d features'%len(params),
        )

        if verbose:
            print('    saving : '+debug_figname)
        fig.savefig(debug_figname)
        plt.close(fig)

        debug_figname = debug_figname.split('.')
        tmp = '.'.join(debug_figname[:-1]) + '-%d.' + debug_figname[-1]

        for ind, region in enumerate(regions):
            fig = data2moi_features_figure(
                points,
                [region],
                rhoc,
                M,
                I,
                baryon_density,
                cs2c2,
                arctan_dlnI_dlnM,
                dlnM_drhoc,
                dlnI_drhoc,
                '%d features'%len(params),
            )

            path = tmp % ind
            if verbose:
                print('    saving : '+path)
            fig.savefig(path)
            plt.close(fig)

    #---

    ### finish formatting data
    params = [[ind]+thing for ind, thing in enumerate(params)] ### include a transition number for reference

    if not len(params):
        params = np.empty((0,len(names)), dtype=float) ### do this so we have a consistent shape

    else:
        params = np.array(params, dtype=float)

    return params, names

#------------------------

def data2moi_features_figure(
        points,
        regions,
        rhoc,
        M,
        I,
        baryon_density,
        cs2c2,
        arctan_dlnI_dlnM,
        dlnM_drhoc,
        dlnI_drhoc,
        title,
    ):
    """make a figure showing the logic encoded in data2moi_features
    """
    fig = plt.figure(figsize=(10, 10))

    axp = plt.subplot(2, 2, 1)
    axc = plt.subplot(2, 2, 2)

    axM1 = plt.subplot(6, 2, 7)
    axM2 = plt.subplot(6, 2, 9)
    axM3 = plt.subplot(6, 2,11)

    axr1 = plt.subplot(6, 2, 8)
    axr2 = plt.subplot(6, 2,10)
    axr3 = plt.subplot(6, 2,12)

    #---

    # plot basic curves
    kwargs = dict(color='k', alpha=0.75)

    axp.plot(rhoc*dlnM_drhoc, rhoc*dlnI_drhoc, **kwargs)
    axc.plot(baryon_density, cs2c2, **kwargs)

    for a, A, y in [
            (axM1, axr1, I),
            (axM2, axr2, arctan_dlnI_dlnM),
            (axM3, axr3, (dlnI_drhoc**2 + dlnM_drhoc**2)**0.5 * rhoc),
        ]:
        a.plot(M, y, **kwargs)
        A.plot(rhoc, y, **kwargs)

    #---

    # plot individual points
    for rho, kwargs in points:
        kwargs = dict(kwargs.items()) ### make a copy so we don't mess with caller

        x = rho*np.interp(rho, rhoc, dlnM_drhoc)
        y = rho*np.interp(rho, rhoc, dlnI_drhoc)

        color = axp.plot(x, y, **kwargs)[0].get_color()
        kwargs['color'] = color

        axc.plot(rho, np.interp(rho, baryon_density, cs2c2), **kwargs)

        for a, A, y in [
                (axM1, axr1, np.interp(rho, rhoc, I)),
                (axM2, axr2, np.interp(rho, rhoc, arctan_dlnI_dlnM)), 
                (axM3, axr3, (x**2 + y**2)**0.5)
            ]:
            a.plot(np.interp(rho, rhoc, M), y, **kwargs)
            A.plot(rho, y, **kwargs)

    #---

    # plot regions

    for rhos, kwargs in regions:
        rhos = sorted(rhos)
        kwargs = dict(kwargs.items()) ### make a copy
        marker = kwargs.pop('marker', '.')

        # add lines
        truth = (rhos[0] <= rhoc) * (rhoc <= rhos[-1])
        lines = [(axp, dlnM_drhoc[truth]*rhoc[truth], dlnI_drhoc[truth]*rhoc[truth])]
        for a, A, y in [
                (axM1, axr1, I),
                (axM2, axr2, arctan_dlnI_dlnM),
                (axM3, axr3, rhoc*(dlnM_drhoc**2 + dlnI_drhoc**2)**0.5)
            ]:
            lines += [(a, M[truth], y[truth]), (A, rhoc[truth], y[truth])]

        truth = (rhos[0] <= baryon_density) * (baryon_density <= rhos[-1])
        lines.append((axc, baryon_density[truth], cs2c2[truth]))

        for a, x, y in lines:
            color = a.plot(x, y, **kwargs)[0].get_color()
            kwargs['color'] = color

        # add point annotations
        
        kwargs['markeredgewidth'] = 2
        kwargs['markeredgecolor'] = color
        kwargs['markerfacecolor'] = 'none'
        kwargs['alpha'] = 0.50

        for points, mark, markersize in [
                (rhos[1:-1], marker, 5),
                ((rhos[0], rhos[-1]), '|', 15),
            ]:
            for rho in points:
                x = rho*np.interp(rho, rhoc, dlnM_drhoc)
                y = rho*np.interp(rho, rhoc, dlnI_drhoc)

                axp.plot(x, y, marker=mark, markersize=markersize, **kwargs)
                axc.plot(rho, np.interp(rho, baryon_density, cs2c2), marker=mark, markersize=markersize, **kwargs)

                for a, A, y in [
                        (axM1, axr1, np.interp(rho, rhoc, I)),
                        (axM2, axr2, np.interp(rho, rhoc, arctan_dlnI_dlnM)),
                        (axM3, axr3, (x**2 + y**2)**0.5)
                    ]:
                    a.plot(np.interp(rho, rhoc, M), y, marker=mark, markersize=markersize, **kwargs)
                    A.plot(rho, y, marker=mark, markersize=markersize, **kwargs)

    #---

    # decorate

    fig.suptitle(title)

#    raise NotImplementedError('add a legend for the different colors/markers used')

    #---

    xlim = axp.get_xlim()
    ylim = axp.get_ylim()

    axp.plot(xlim, [0]*2, color='k', alpha=0.1)
    axp.plot([0]*2, ylim, color='k', alpha=0.1)

    if xlim[0] < 0:
        axp.fill_between([xlim[0], 0], [ylim[0]]*2, [ylim[1]]*2, color='k', alpha=0.05)

    axp.set_xlim(xlim)
    axp.set_ylim(ylim)

    axp.set_xlabel('$d\ln M/d\ln \\rho_c$')
    axp.set_ylabel('$d\ln I/d\ln \\rho_c$')

    axp.xaxis.tick_top()
    axp.xaxis.set_label_position('top')

    #---

    axc.set_ylabel('$c_s^2/c^2$')
    axc.set_xlabel('$\\rho$')
    axc.xaxis.tick_top()
    axc.xaxis.set_label_position('top')

    axc.yaxis.tick_right()
    axc.yaxis.set_label_position('right')

    axc.set_xscale('log')
    axc.set_yscale('linear')

    axc.set_ylim(ymin=-0.05, ymax=1.05)

    xlim = axc.get_xlim()
    if xlim[0] < 9.0e13: ### 0.1*rho_nuc
        axc.set_xlim(xmin=9.0e13)
#        axc.set_ylim(ymin=np.min(cs2c2[baryon_density>=2.8e13]))

    if xlim[1] > 1.1*rhoc[-1]:
        axc.set_xlim(xmax=1.1*rhoc[-1])

    axc.grid(True, which='both')

    #---

    axM1.set_ylabel('$I$')
    axr1.set_ylabel(axM1.get_ylabel())
    axr1.yaxis.tick_right()
    axr1.yaxis.set_label_position('right')

    plt.setp(axM1.get_xticklabels(), visible=False)
    plt.setp(axr1.get_xticklabels(), visible=False)

    axr1.set_xscale('log')
    axr1.set_xlim(axc.get_xlim())

    axM1.grid(True, which='both')
    axr1.grid(True, which='both')

    #---

    axM2.set_ylabel('$\\tan^{-1}\left(\\frac{d\ln I}{d\ln M}\\right)$')
    axr2.set_ylabel(axM2.get_ylabel())
    axr2.yaxis.tick_right()
    axr2.yaxis.set_label_position('right')

    plt.setp(axM2.get_xticklabels(), visible=False)
    plt.setp(axr2.get_xticklabels(), visible=False)

    axM2.set_xscale(axM1.get_xscale())
    axM2.set_xlim(axM1.get_xlim())

    axr2.set_xscale(axr1.get_xscale())
    axr2.set_xlim(axr1.get_xlim())

    axM2.grid(True, which='both')
    axr2.grid(True, which='both')

    for ax2 in [axM2, axr2]:
        ylim = ax2.get_ylim()
        if ylim[0] < -0.5*np.pi:
            ax2.fill_between(ax2.get_xlim(), [ylim[0]]*2, [-0.5*np.pi]*2, color='k', alpha=0.05)
        if ylim[1] > +0.5*np.pi:
            ax2.fill_between(ax2.get_xlim(), [ylim[1]]*2, [+0.5*np.pi]*2, color='k', alpha=0.05)
        ax2.set_ylim(ylim)

    #---

#    axM3.set_ylabel('$\sqrt{(d\ln I/d\ln \\rho_c)^2 + (d\ln M/d\ln\\rho_c)^2}$')
    axM3.set_ylabel('$\sqrt{(d\ln I)^2 + (d\ln M)^2}$')
    axr3.set_ylabel(axM3.get_ylabel())
    axr3.yaxis.tick_right()
    axr3.yaxis.set_label_position('right')

    axM3.set_yscale('log')
    axr3.set_yscale('log')

    axM3.set_xlabel('$M$')
    axr3.set_xlabel(axc.get_xlabel())

    axM3.set_xscale(axM2.get_xscale())
    axM3.set_xlim(axM2.get_xlim())

    axr3.set_xscale(axr2.get_xscale())
    axr3.set_xlim(axr2.get_xlim())

    axM3.grid(True, which='both')
    axr3.grid(True, which='both')

    #---

    plt.subplots_adjust(
        left=0.10,
        right=0.90,
        top=0.92,
        bottom=0.05,
        hspace=0.08,
        wspace=0.03,
    )

    return fig

#------------------------

def process2moi_features(
        data,
        eos_template,
        eos_num_per_dir,
        mac_template,
        macro_num_per_dir,
        summary_template,
        eos_rho,
        eos_cs2c2,
        macro_rhoc,
        macro_mass,
        macro_moi,
        output_eos_columns=DEFAULT_EOS_COLUMNS,
        output_macro_columns=DEFAULT_MACRO_COLUMNS,
        flatten_thr=DEFAULT_FLATTEN_THR,
        smoothing_width=DEFAULT_SMOOTHING_WIDTH,
        diff_thr=DEFAULT_DIFF_THR,
        cs2c2_drop_ratio=DEFAULT_CS2C2_DROP_RATIO,
        verbose=False,
        debug=False,
    ):
    """extract the branches for each EoS specified
    """
    N = len(data)
    for ind, eos in enumerate(data):

        ### construct paths
        # where we're going to read in data
        eos_path = eos_template%{'moddraw':eos//eos_num_per_dir, 'draw':eos}

        tmp = {'moddraw':eos//macro_num_per_dir, 'draw':eos}
        mac_path = mac_template%tmp

        # where we're going to write data
        sum_path = summary_template%tmp

        if verbose:
            print('    %d/%d'%(ind+1, N))
            print('    loading macro: %s'%mac_path)
        mac_data, mac_cols = io.load(mac_path, [macro_rhoc, macro_mass, macro_moi]+output_macro_columns) ### NOTE: we load all columns because we're going to re-write them all into subdir as separate branches

        if verbose:
            print('    loading eos: %s'%eos_path)
        eos_data, eos_cols = io.load(eos_path, [eos_rho, eos_cs2c2]+output_eos_columns) ### NOTE: this guarantees that eos_rho is the first column!
        baryon_density = eos_data[:,eos_cols.index(eos_rho)] ### separate this for convenience
        cs2c2 = eos_data[:,eos_cols.index(eos_cs2c2)]

        # use macro data to identify separate stable branches
        # NOTE: we expect this to be ordered montonically in rhoc
        M = mac_data[:,mac_cols.index(macro_mass)]
        I = mac_data[:,mac_cols.index(macro_moi)]
        rhoc = mac_data[:,mac_cols.index(macro_rhoc)]

        if debug:
            debug_figname = eos_path[:-4] + '-moi_features.png'
        else:
            debug_figname = None

        params, names = data2moi_features(
            rhoc,
            M,
            I,
            baryon_density,
            cs2c2,
            mac_data,
            mac_cols,
            eos_data,
            eos_cols,
            flatten_thr=flatten_thr,
            smoothing_width=smoothing_width,
            diff_thr=diff_thr,
            cs2c2_drop_ratio=cs2c2_drop_ratio,
            verbose=verbose,
            debug_figname=debug_figname,
        )

        if verbose:
            print('    writing summary of %d identified moi-features into: %s'%(len(params), sum_path))

        newdir = os.path.dirname(sum_path)
        if not os.path.exists(newdir):
            try:
                os.makedirs(newdir)
            except OSError:
                pass ### directory already exists

        io.write(sum_path, params, names)



# andal's functions here 

from scipy.ndimage import gaussian_filter1d

def tanh_k_transform(arctan_data, central_pressure, mass, sw=3.5):
    """compute the tanh K transform for additional features based on the transformed mass-radius relationship
    """

    dm_dpc = utils.num_dfdx(np.log(central_pressure), (mass)) 
    dm_dpc_second = utils.num_dfdx(np.log(central_pressure), dm_dpc) 
    darctan_D_dpc = utils.num_dfdx(np.log(central_pressure), arctan_data)
    darctan_D_dpc_second = utils.num_dfdx(np.log(central_pressure), darctan_D_dpc) 
    raw_k_r = (dm_dpc * darctan_D_dpc_second -  darctan_D_dpc * dm_dpc_second) / ( dm_dpc**2 + darctan_D_dpc**2)**1.5

    smoothed_k_r = gaussian_filter1d(raw_k_r, sigma=sw) # gaussian smoothing to filter out noise, ~ 3.5 seems to work well 

    return np.tanh(smoothed_k_r) 


def modify_ticks(plot):
    plot.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    plot.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))

    for label in plot.get_xticklabels() + plot.get_yticklabels():
        label.set_fontsize(30)
        label.set_fontweight('extra bold')
        label.set_color('black')

    plot.xaxis.label.set_size(32)
    plot.yaxis.label.set_size(32)
    plot.title.set_size(18)
    plot.title.set_fontweight(900)       

    plot.xaxis.label.set_fontweight('extra bold')
    plot.yaxis.label.set_fontweight(1000)

    plot.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    plot.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    plot.tick_params(axis='both', which='major', width=2.5, length=10, color='black', direction='in', pad=10)
    plot.tick_params(top=True, bottom=True, left=True, right=True)

    for spine in plot.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(2.5)


DEFAULT_K_FEATURE_NAME = "feature"
DEFAULT_ARCTAN_DIFF_THRESHOLD= 0.05
DEFAULT_TANHK_DIFF_THRESHOLD = 0.05
DEFAULT_CS2C2_DROP_RATIO = 1.11
MAX_K_TEMPLATE = "max_tanh_k_%s"
MIN_K_TEMPLATE = "min_tanh_k_%s"   

def data2_tanh_features(
          rhoc,
          mass, 
          radius, 
          baryon_density, 
          eos_cs2c2, 
          macro_data,
          macro_cols, 
          eos_data, 
          eos_cols, 
          central_pressure,
          verbose=True,
          flatten_thr=DEFAULT_FLATTEN_THR,
          smoothing_width=DEFAULT_SMOOTHING_WIDTH,
          diff_thr=DEFAULT_ARCTAN_DIFF_THRESHOLD,
          cs2c2_drop_ratio=DEFAULT_CS2C2_DROP_RATIO,
          tanh_k_threshold=DEFAULT_TANHK_DIFF_THRESHOLD, 
          tanh_k_sigma=3,
          cs_drop_threshold=DEFAULT_CS2C2_DROP_RATIO, 
          diff_k_threshold=DEFAULT_TANHK_DIFF_THRESHOLD,
          diff_arctan_threshold=DEFAULT_ARCTAN_DIFF_THRESHOLD
    ):

    names = [DEFAULT_K_FEATURE_NAME]

    for tmp in [
            MIN_CS2C2_TEMPLATE, 
            MAX_CS2C2_TEMPLATE,
            RMAX_CS2C2_TEMPLATE,
            LOW_CS2C2_TEMPLATE,
            HGH_CS2C2_TEMPLATE,
            MIN_K_TEMPLATE,
            MAX_K_TEMPLATE, 
        ]:
        names.append(tmp % 'tanh_k')
        names.append(tmp % 'arctan_dlnR_dlnM')
        names += [tmp%col for col in macro_cols]
        names += [tmp%col for col in eos_cols]
    
    params= []

    arctan_dlnR_dlnM, (spurious, dlnM_drhoc, dlnI_drhoc) = arctan_transform(rhoc, mass, radius, smoothing_width=0.01) 

    tanh_k = tanh_k_transform(arctan_dlnR_dlnM, rhoc, mass, sw = tanh_k_sigma) 
    maxima_in_tanh_k = find_inclusive_maxima(tanh_k)


    ends = list([index for index in maxima_in_tanh_k if tanh_k[index] > 0][::-1]) # new indicator variable
    ends += list(find_minima(arctan_dlnR_dlnM)[::-1]) # add features from old indicator variable as well 

    while len(ends):
        end = ends[-1]
        if np.any(dlnM_drhoc[:end] > 0): ### something is stable before this
             break
        ends = ends[:-1] ### truncate this guy

    # discard any local minima that are in the final unstable branch
    while len(ends):
        end = ends.pop(0)
        if np.any(dlnM_drhoc[end:] > 0): ### not part of the 'final unstable branch'
            ends.insert(0, end)
            break

    if ends: ### we have something to do

        ### if the lowest rhoc available is a minimum, we discard it since we can't tell if it is a local minimum or not

        if ends[-1] == 0:                                                  # ANDAL: this is basically that if the min in arctan(D) is due to the first point, so ignore that. 
            ends = ends[:-1]

        if ends[0] == len(rhoc)-1: ### same thing with the end point
            ends = ends[1:]
        
            ### local minima in sound speed
        min_cs2c2 = np.array(find_inclusive_minima(eos_cs2c2))
        min_cs2c2_baryon_density = baryon_density[min_cs2c2]

        ### local maxima in sound speed
        max_cs2c2 = np.array(find_inclusive_maxima(eos_cs2c2))
        max_cs2c2_baryon_density = baryon_density[max_cs2c2]

        rmax_cs2c2 = find_running_maxima(eos_cs2c2)

        if not rmax_cs2c2: ### no qualifying maxima
                
                ends = [] ### this will make us skip all the ends because we can't match them to starting points

        else:
                rmax_cs2c2 = np.array(rmax_cs2c2)
                rmax_cs2c2_baryon_density = baryon_density[rmax_cs2c2]

        Neos = len(eos_cols)
        Nmac = len(macro_cols)

        mtov = np.argmax(mass)

        for end in ends:
                r = rhoc[end]

                mtov_baryon_density = rhoc[mtov]

                if r >= mtov_baryon_density:
                    continue


                if verbose:
                    print('        processing end=%d/%d at rhoc=%.6e'%(end, len(rhoc), r))

                #---

                ### min sound speed preceeding "end"
                try:
                    ind_min_cs2c2 = find_preceeding(r, min_cs2c2_baryon_density, min_cs2c2)
                except RuntimeError:
                    if verbose:
                        print('            WARNING! could not find local minimum in cs2c2 preceeding local minimum in arctan(dlogI/dlogM)')

                    continue

                min_r = baryon_density[ind_min_cs2c2]
                min_cs2c2_tanh_k = np.interp(min_r, rhoc, tanh_k)

                #---

                ### max sound speed preceeding minimum sound speed
                try:
                    ind_max_cs2c2 = find_preceeding(min_r, max_cs2c2_baryon_density, max_cs2c2)
                except RuntimeError:
                    if verbose:
                        print('            WARNING! could not find local maximum in cs2c2 preceeding local minimum in cs2c2')

                    continue

                max_r = baryon_density[ind_max_cs2c2]

                try:
                    ind_rmax_cs2c2 = find_preceeding(max_r, rmax_cs2c2_baryon_density, rmax_cs2c2)
                except RuntimeError:
                    if verbose:
                        print('            WARNING! could not find running maximum in cs2c2 preceeding local minimum in cs2c2')

                    continue

                rmax_r = baryon_density[ind_rmax_cs2c2]
                rmax_cs2c2_tanh_k = np.interp(rmax_r, rhoc, tanh_k)  #ANDAL: replaced arctan with tanh(K) here 
                
                selected = (rmax_r <=rhoc) * (rhoc <= r)
                max_tanh_k = rhoc[selected][np.argmax(tanh_k[selected])]
                ind_max_tanh_k = np.arange(len(rhoc))[rhoc==max_tanh_k][0]

                max_arctan_r = rhoc[selected][np.argmax(arctan_dlnR_dlnM[selected])]
                ind_max_arctan= np.arange(len(rhoc))[rhoc==max_arctan_r][0]

                selected = (rmax_r <= baryon_density) * (baryon_density<=r)
                low_cs2c2_r = baryon_density[selected][np.argmin(eos_cs2c2[selected])]
                ind_low_cs2c2 = np.arange(len(eos_cs2c2))[low_cs2c2_r==baryon_density][0]

                hgh_cs2c2_r = baryon_density[selected][np.argmax(eos_cs2c2[selected])]
                ind_hgh_cs2c2 = np.arange(len(eos_cs2c2))[hgh_cs2c2_r==baryon_density][0]


                ####

                diff_k = rmax_cs2c2_tanh_k - tanh_k[end]  # this should be difference between running max and end in tanhk to remove small "wiggles" in tanhk 
                diff_arctan = arctan_dlnR_dlnM[ind_max_arctan] - arctan_dlnR_dlnM[end]

                # sound speed must drop by a certain fraction
                drop_ratio = eos_cs2c2[ind_rmax_cs2c2] / eos_cs2c2[ind_low_cs2c2]
                cs2c2_drop_ratio = cs_drop_threshold #ANDAL: why is the default 0..?
                DIFF_K_THRESH = diff_k_threshold
                DIFF_ARCTAN_RATIO = diff_arctan_threshold

                while ((drop_ratio < cs2c2_drop_ratio) or (diff_arctan < DIFF_ARCTAN_RATIO) or (np.abs(diff_k) < DIFF_K_THRESH)): # post-hoc removal of phase transitions 

                    # first, find preceeding local max_cs2c2
                    try:
                        truth = max_cs2c2_baryon_density < rmax_r ### move to a local max that is before the current running max
                        jnd_max_cs2c2 = find_preceeding(rmax_r, max_cs2c2_baryon_density[truth], max_cs2c2[truth])

                    except RuntimeError:
                        if verbose:
                            print('            WARNING! could not find local maximum in cs2c2 preceeding local minimum in cs2c2')
                        break

                    # then update rmax_cs2c2 to be before that
                    try:
                        ind_rmax_cs2c2 = find_preceeding(baryon_density[jnd_max_cs2c2], rmax_cs2c2_baryon_density, rmax_cs2c2)
                    except RuntimeError:
                        if verbose:
                            print('            WARNING! could not find running maximum in cs2c2 preceeding local minimum in cs2c2')

                        break

                    rmax_r = baryon_density[ind_rmax_cs2c2]
                    rmax_cs2c2_tanh_k = np.interp(rmax_r, rhoc, tanh_k)

                    ### now find the maximum arctan(...) between rmax_r and r
                    selected = (rmax_r<=rhoc) * (rhoc<=r)
                    max_tanh_k = rhoc[selected][np.argmax(tanh_k[selected])]
                    ind_max_tanh_k = np.arange(len(rhoc))[rhoc==max_tanh_k][0]

                    max_arctan_r = rhoc[selected][np.argmax(arctan_dlnR_dlnM[selected])]
                    ind_max_arctan= np.arange(len(rhoc))[rhoc==max_arctan_r][0]

                    selected = (rmax_r<=baryon_density) * (baryon_density<=r)
                    low_cs2c2_r = baryon_density[selected][np.argmin(eos_cs2c2[selected])]
                    ind_low_cs2c2 = np.arange(len(eos_cs2c2))[low_cs2c2_r==baryon_density][0]

                    hgh_cs2c2_r = baryon_density[selected][np.argmin(eos_cs2c2[selected])]
                    ind_hgh_cs2c2 = np.arange(len(eos_cs2c2))[hgh_cs2c2_r==baryon_density][0]

                    
                    # update conditionals
                    diff_k = rmax_cs2c2_tanh_k - tanh_k[end]
                    diff_arctan = arctan_dlnR_dlnM[ind_max_arctan] - arctan_dlnR_dlnM[end]
                    drop_ratio = eos_cs2c2[ind_rmax_cs2c2] / eos_cs2c2[ind_low_cs2c2]

                if (drop_ratio < cs2c2_drop_ratio) or (diff_arctan < DIFF_ARCTAN_RATIO) or (np.abs(diff_k) < DIFF_K_THRESH): # we exited from the break statement
                    continue

                datum = []

                # add local min in cs2c2
                datum.append(np.interp(min_r, rhoc, tanh_k))
                datum.append(np.interp(min_r, rhoc, arctan_dlnR_dlnM))
                datum += [np.interp(min_r, rhoc, macro_data[:,i]) for i in range(Nmac)]
                datum += list(eos_data[ind_min_cs2c2])

                # add local max in cs2c2
                datum.append(np.interp(max_r, rhoc, tanh_k))
                datum.append(np.interp(max_r, rhoc, arctan_dlnR_dlnM))
                datum += [np.interp(max_r, rhoc, macro_data[:,i]) for i in range(Nmac)]
                datum += list(eos_data[ind_max_cs2c2])

                # add running max in cs2c2
                datum.append(np.interp(rmax_r, rhoc, tanh_k))
                datum.append(np.interp(rmax_r, rhoc, arctan_dlnR_dlnM))
                datum += [np.interp(rmax_r, rhoc, macro_data[:,i]) for i in range(Nmac)]
                datum += list(eos_data[ind_rmax_cs2c2])

                # add lowest cs2c2 within feature
                datum.append(np.interp(low_cs2c2_r, rhoc, tanh_k))
                datum.append(np.interp(low_cs2c2_r, rhoc, arctan_dlnR_dlnM))
                datum += [np.interp(low_cs2c2_r, rhoc, macro_data[:,i]) for i in range(Nmac)]
                datum += list(eos_data[ind_low_cs2c2])

                # add highest cs2c2 within feature
                datum.append(np.interp(hgh_cs2c2_r, rhoc, tanh_k))
                datum.append(np.interp(hgh_cs2c2_r, rhoc, arctan_dlnR_dlnM))
                datum += [np.interp(hgh_cs2c2_r, rhoc, macro_data[:,i]) for i in range(Nmac)]
                datum += list(eos_data[ind_hgh_cs2c2])

                # add parameters at max tanhk
                datum.append(tanh_k[ind_max_tanh_k])
                datum.append(arctan_dlnR_dlnM[ind_max_tanh_k])
                datum += list(macro_data[ind_max_tanh_k])
                datum += [np.interp(max_tanh_k, baryon_density, eos_data[:,i]) for i in range(Neos)]

                # add parameters at "end" (min tanhk)
                datum.append(tanh_k[end])
                datum.append(arctan_dlnR_dlnM[end])
                datum += list(macro_data[end])
                datum += [np.interp(r, baryon_density, eos_data[:,i]) for i in range(Neos)]

                params.append(datum)

        
        ### finish formatting data
        params = [[ind]+thing for ind, thing in enumerate(params)] ### include a transition number for reference

        if not len(params):
            params = np.empty((0,len(names)), dtype=float) ### do this so we have a consistent shape

        else:
            params = np.array(params, dtype=float)

    return (params, names, arctan_dlnR_dlnM, tanh_k)


def process2all_features(
        data,                   
        eos_template,
        eos_num_per_dir,
        mac_template,
        macro_num_per_dir,
        summary_template,
        eos_rho,
        eos_cs2c2,
        macro_rhoc,
        macro_mass,
        macro_radius,
        output_eos_columns=DEFAULT_EOS_COLUMNS,
        output_macro_cols=DEFAULT_MACRO_COLUMNS,
        flatten_thr=DEFAULT_FLATTEN_THR,
        smoothing_width=DEFAULT_SMOOTHING_WIDTH,
        diff_thr=DEFAULT_DIFF_THR,
        cs2c2_drop_ratio=DEFAULT_CS2C2_DROP_RATIO,
        diff_arctan_threshold=DEFAULT_ARCTAN_DIFF_THRESHOLD,
        diff_k_threshold=DEFAULT_TANHK_DIFF_THRESHOLD,
        verbose=False,
        debug=False,
    ):
    """extract the branches for each EoS specified
    """
    N = len(data)                                                              
    for ind, eos in enumerate(data): # so it does have some eos.. 

        ### construct paths
        # where we're going to read in data
        eos_path = eos_template%{'moddraw':eos//eos_num_per_dir, 'draw':eos}   

        tmp = {'moddraw':eos//macro_num_per_dir, 'draw':eos}
        mac_path = mac_template%tmp

        # where we're going to write data
        sum_path = summary_template%tmp

        if verbose:
            print('    %d/%d'%(ind+1, N))
            print('    loading macro: %s'%mac_path)
        mac_data, mac_cols = io.load(mac_path, [macro_rhoc, macro_mass, macro_radius]+output_macro_cols) ### NOTE: we load all columns because we're going to re-write them all into subdir as separate branches

        if verbose:
            print('    loading eos: %s'%eos_path)
        eos_data, eos_cols = io.load(eos_path, [eos_rho, eos_cs2c2]+output_eos_columns) ### NOTE: this guarantees that eos_rho is the first column!
        baryon_density = eos_data[:,eos_cols.index(eos_rho)] ### separate this for convenience
        eos_cs2c2 = eos_data[:,eos_cols.index(eos_cs2c2)] 

        # use macro data to identify separate stable branches
        # NOTE: we expect this to be ordered montonically in rhoc
        M = mac_data[:,mac_cols.index(macro_mass)]
        radius = mac_data[:,mac_cols.index(macro_radius)] 
        rhoc = mac_data[:,mac_cols.index(macro_rhoc)]
        central_pressure = mac_data[:,mac_cols.index("central_pressurec2")]

        if debug:
            debug_figname = eos_path[:-4] + '-all_features.png'
        else:
            debug_figname = None

        params, names, arctan_dlnR_dlnM, tanh_k = data2_tanh_features(
            rhoc, 
            mass=M, 
            radius=radius, 
            eos_cs2c2=eos_cs2c2, 
            baryon_density=baryon_density,
            central_pressure=central_pressure,
            macro_data= mac_data, 
            eos_data=eos_data, 
            eos_cols = eos_cols, 
            macro_cols=mac_cols,
            diff_arctan_threshold=diff_arctan_threshold, 
            cs_drop_threshold=cs2c2_drop_ratio,
            diff_k_threshold=diff_k_threshold)
        
        if verbose:
            print('    writing summary of %d identified all-features into: %s'%(len(params), sum_path))

        newdir = os.path.dirname(sum_path)
        if not os.path.exists(newdir):
            try:
                os.makedirs(newdir)
            except OSError:
                pass ### directory already exists

        io.write(sum_path, params, names)


