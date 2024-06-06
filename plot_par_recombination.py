#!/usr/bin/env python

# Scrapes recombination rate files
# Averages on a per-base level across PAR1
# Plots results and saves as .png and .svg
# Also calculates and plots cumulative
# recombination frequency on an averaged
# and per donor basis.

# Packages
import argparse
import os
import fnmatch
import re
from math import isclose
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FuncFormatter
from matplotlib import font_manager as fm

# Global variables
telomeric_boundary = 251087
proposed_PAB = 2281479
established_PAB = 2781479

hinch_proximal_marker = 2352082
hinch_distal_marker = 2781479

# Utility Functions

def find_files(pattern, path):
    '''Returns list of all recombination freq files'''
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return result

def custom_key(item):
    '''Function reorders donors lexicographically'''
    return int(re.search(r'\d+', item).group())

# Data Manipulation Functions

def get_recomb_array(df):
    '''Returns 1D numpy array with recombination frequencies'''
    start = df['Start'].min() # 1-indexed
    end = df['End'].max() # 1-indexed
    # create 1D array with np.nan
    recomb_array = np.full((end-start), fill_value=np.NaN)
    # Fill in array by row
    for row in df.itertuples(index=True, name='Pandas'):
        # get start location on array
        for i in (range(row.Start-start, row.End-start)):
            # Leaving 0 recombination rates np.nan since not observed
            recomb_freq = getattr(row, 'RF')
            if (recomb_freq != 0):
                recomb_array[i] = recomb_freq
    return recomb_array

def get_cum_AUC_proximal_from_telomere(df):
    '''
    Given Recombination Frequency Array return
    cumulative recombination frequency AUC
    with using the telomere as the start point.
    Inputs:
    df := Pandas DataFrame, recombination frequencies at each bp
    '''
    dx = 1/10**6 # Get interval
    # Define holder array for AUC
    AUC_array_cumulative = []
    PAB_AUC_cumulative = []
    # Cycle through donors
    for i, donor in enumerate(df.columns):
        # Calculating AUC across entire PAR1
        genomic_coordinates = df.index/10**6
        mask = ~np.isnan(df[donor])
        genetic_distance = np.full(len(genomic_coordinates),
                                   fill_value=np.nan)
        genetic_distance[0] = 0 # Initialize count to 0
        # Iterate through recombination frequencies and calculate running AUC
        for j in range(1, len(genetic_distance)):
            if mask[j] == True:
                # Add the area of the rectangle formed by current RF
                # and the previous RF
                genetic_distance[j] = genetic_distance[j - 1] + (df[donor][j] * dx)
            else:
                genetic_distance[j] = genetic_distance[j - 1]
        AUC_array_cumulative.append(genetic_distance)

        # Calculating AUC only for bases between PABs
        genomic_coordinates = (df.index/10**6)[proposed_PAB:established_PAB]
        PAB_interval_RF = df[donor][proposed_PAB:established_PAB]
        mask = ~np.isnan(PAB_interval_RF)
        genetic_distance = np.full(len(genomic_coordinates),
                                   fill_value=np.nan)
        genetic_distance[0] = 0 # Initialize count to 0
        # Iterate through recombination frequencies and calculate running AUC
        for h, j in enumerate(range(np.min(PAB_interval_RF.index), np.max(PAB_interval_RF.index)+1)):
            if mask[j] == True and h == 0:
                genetic_distance[h] = PAB_interval_RF[j] * dx
            elif mask[j] == True and h != 0:
                # Add the area of the rectangle formed by current RF
                # and the previous RF
                genetic_distance[h] = genetic_distance[h - 1] + (PAB_interval_RF[j] * dx)
            else:
                genetic_distance[h] = genetic_distance[h - 1]
        PAB_AUC_cumulative.append(genetic_distance)
    # Convert list of arrays to dataframe
    AUC_array_cumulative = pd.DataFrame(list(map(np.ravel, AUC_array_cumulative)))
    AUC_array_cumulative = AUC_array_cumulative.T
    AUC_array_cumulative.columns = df.columns
    PAB_AUC_cumulative = pd.DataFrame(list(map(np.ravel, PAB_AUC_cumulative)))
    PAB_AUC_cumulative = PAB_AUC_cumulative.T
    PAB_AUC_cumulative.columns = df.columns
    return AUC_array_cumulative, PAB_AUC_cumulative

def get_cum_AUC_distal_from_PAR(df):
    '''
    Given Recombination Frequency Array return
    cumulative recombination frequency AUC
    with using the established PAR1 boundary
    as the start point.
    Inputs:
    df := Pandas DataFrame, recombination frequencies at each bp
    '''
    dx = 1/10**6 # Get interval
    # Define holder array for AUC
    AUC_array_cumulative = []
    PAB_AUC_cumulative = []
    # Cycle through donors
    for i, donor in enumerate(df.columns):
        # Calculating AUC across entire PAR1
        genomic_coordinates = df.index/10**6
        mask = ~np.isnan(df[donor])
        genetic_distance = np.full(len(genomic_coordinates),
                                   fill_value=np.nan)
        genetic_distance[-1] = 0 # Initialize count to 0
        # Iterate through recombination frequencies and calculate running AUC
        for j in reversed(range(0, len(genetic_distance)-1)):
            if mask[j] == True:
                # Add the area of the rectangle formed by current RF
                # and the previous RF
                genetic_distance[j] = genetic_distance[j + 1] + (df[donor][j] * dx)
            else:
                genetic_distance[j] = genetic_distance[j + 1]
        AUC_array_cumulative.append(genetic_distance)

        # Calculating AUC only for bases between PABs
        genomic_coordinates = (df.index/10**6)[proposed_PAB:established_PAB]
        PAB_interval_RF = df[donor][proposed_PAB:established_PAB]
        mask = ~np.isnan(PAB_interval_RF)
        genetic_distance = np.full(len(genomic_coordinates),
                                   fill_value=np.nan)
        genetic_distance[0] = 0 # Initialize count to 0
        # Iterate through recombination frequencies and calculate running AUC
        for h, j in enumerate(reversed(range(np.min(PAB_interval_RF.index), np.max(PAB_interval_RF.index)+1))):
            if mask[j] == True and h == 0:
                genetic_distance[h] = PAB_interval_RF[j] * dx
            elif mask[j] == True and h != 0:
                # Add the area of the rectangle formed by current RF
                # and the previous RF
                genetic_distance[h] = genetic_distance[h - 1] + (PAB_interval_RF[j] * dx)
            elif h != 0:
                genetic_distance[h] = genetic_distance[h - 1]
        PAB_AUC_cumulative.append(genetic_distance)
    # Convert list of arrays to dataframe
    AUC_array_cumulative = pd.DataFrame(list(map(np.ravel, AUC_array_cumulative)))
    AUC_array_cumulative = AUC_array_cumulative.T
    AUC_array_cumulative.columns = df.columns
    PAB_AUC_cumulative = pd.DataFrame(list(map(np.ravel, PAB_AUC_cumulative)))
    PAB_AUC_cumulative = PAB_AUC_cumulative.T
    PAB_AUC_cumulative.columns = df.columns
    PAB_AUC_cumulative = PAB_AUC_cumulative[::-1].reset_index(drop=True) # Flip index
    return AUC_array_cumulative, PAB_AUC_cumulative

def check_AUC(df_proximal, df_distal, portion_of_par):
    '''
    Makes sure AUC calculations are same going from proximal
    to distal direction.
    Inputs:
    df_proximal := AUC dataframe calculated going from telomere
    df_distal := AUC dataframe calculated going from PAR1 boundary
    portion_of_par := string, what part of PAR1 {entirety, contested boundary}
    '''
    print('Checking {} of the PAR:'.format(portion_of_par))
    max_AUC_proximal = [x for x in df_proximal.max()]
    max_AUC_distal = [x for x in df_distal.max()]
    truth_array = []
    for i in range(len(max_AUC_proximal)):
        truth_array.append(isclose(max_AUC_proximal[i],
                                   max_AUC_distal[i],
                                   abs_tol=0.001,
                                   ))
    if np.any(truth_array) == False:
        false_indices = np.where(truth_array == False)[0]
        print('Issues at the following donors: {}'.format(df_proximal[false_indices].to_list()))
    else:
        print('Directionality of AUC calculations is consistent across all donors.')
    return

# Plotting Functions

def format_ticks_for_donor_plots(x, pos):
    '''Only shows specific tick values for donor plots'''
    if x == 0 or x == 25 or x == 50:
        return int(x)
    else:
        return ''

def plot_averaged_rf(df, path):
    '''
    Plots average recombination rate across PAR1
    Inputs:
    df := Pandas DataFrame; RF at every basepair
    path := str; save directory
    '''
    # Get limits normalized in megabases (Mb)
    start, end = 1/10**6, df.shape[0]/10**6
    fig, ax = plt.subplots(figsize=(10, 5), dpi=300)
    # Removing problematic recombination rates (>100%)
    df = df.copy()
    df = df.map(lambda x: np.nan if x >= 100 else x)
    # Plotting average
    ax.plot((df.index/10**6)[:proposed_PAB], df.mean(axis=1)[:proposed_PAB],
            '-', c='#B3D4BE', label='Average (PAR1)',
            alpha=1.0,
            )
    ax.plot((df.index/10**6)[proposed_PAB:], df.mean(axis=1)[proposed_PAB:],
            '-', c='#74B088', label='Average (Contested PAB Interval)',
            alpha=1.0,
            )
    # Plotting proposed PAR boundary (PAB)
    plt.axvline(x=proposed_PAB/10**6, ymax=0.95,
                c='black', linestyle=(0, (1, 3, 6, 3, 1, 3)),
                linewidth=1.5,
                label='Proposed PAB',
                )
    # Plotting established PAR boundary (PAB)
    plt.axvline(x=established_PAB/10**6, ymax=0.95,
                c='black', linestyle=(0, (9, 3)),
                linewidth=1.5,
                label='Established PAB',
                )
    # Set graph limits
    ax.set_xlim(left=start, right=end+0.002)
    ax.set_ylim(bottom=0, top=50)
    # Annotate PAB
    ax.annotate('Proposed\nPAB', ha='center',
                xy=(proposed_PAB/10**6, 40), xytext=(proposed_PAB/10**6, 48),
                )
    ax.annotate('Established\nPAB', ha='center',
                xy=(established_PAB/10**6, 40), xytext=(established_PAB/10**6, 48),
                )
    # Labels
    plt.title('Averaged Observed Recombination Rate on PAR1')
    ax.set_ylabel('Rate (cM/Mb)')
    ax.set_xlabel('Genomic Coordinates (Mb)')
    # Clean up graph
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_major_locator(MultipleLocator(10))
    # Changing X-axis to match other figures
    ax.spines['bottom'].set_color('#B3D3BE')
    ax.spines['bottom'].set_linewidth(5)
    ax.tick_params(axis='x', color='white', width=1,
                    direction='inout', length=5,
                    )
    ax.spines['bottom'].set_position(('outward', 5))
    ax.spines['bottom'].set_zorder(0)
    # Plotting
    fig.tight_layout()
    plt.savefig(path+'average_PAR_recomb.svg',
                transparent=False,)
    plt.savefig(path+'average_PAR_recomb.png',
            transparent=False,)
    return

def plot_donor_rf(df, path, ncols=2, figsize=(8, 10)):
    '''
    Plots individual donor level recombination frequency
    at each base pair in the PAR1 region.
    Inputs:
    df := pandas DataFrame; includes RF at each bp position
    path := string; save directory
    ncols := int; number of columns to use in output graph
    figsize := tuple of ints; (width, height) in inches
    '''
    # Get limits noramlized in megabases (Mb)
    start, end = 1/10**6, df.shape[0]/10**6
    # Create figure
    fig = plt.figure(figsize=figsize, dpi=300)
    # Define number of rows / columns for subplot grid
    nrows = df.shape[1] // ncols + (df.shape[1] % ncols > 0)
    for i, donor in enumerate(df.columns):
        ax = plt.subplot(nrows, ncols, (i+1)) # new graph
        # Plotting
        ax.plot((df.index/10**6)[:proposed_PAB], df[donor][:proposed_PAB],
                '-', c='#B3D4BE', label=donor,
                alpha=1.0,
                )
        ax.plot((df.index/10**6)[proposed_PAB:], df[donor][proposed_PAB:],
                '-', c='#74B088', label=donor,
                alpha=1.0,
                )
        # Plotting proposed PAR boundary (PAB)
        plt.axvline(x=proposed_PAB/10**6, ymax=1.00,
                    c='black', linestyle=(0, (1, 3, 6, 3, 1, 3)),
                    linewidth=1.5,
                    label='Proposed PAB',
                    )
        # Plotting established PAR boundary (PAB)
        plt.axvline(x=established_PAB/10**6, ymax=1.00,
                    c='black', linestyle=(0, (9, 3)),
                    linewidth=1.5,
                    label='Established PAB',
                    )
        # Set graph limits
        ax.set_xlim(left=start, right=end+0.02)
        ax.set_ylim(bottom=0, top=50)
        # Labels
        ax.set_title((' ' + donor), loc='left', y=0.90, fontsize=12)
        # Clean up graph
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        # Set major tick locator for y-axis
        ax.yaxis.set_major_locator(MultipleLocator(5))
        ax.yaxis.set_major_formatter(FuncFormatter(format_ticks_for_donor_plots))
        # Set tick locators for x-axis
        ax.xaxis.set_major_locator(MultipleLocator(1.0))
        ax.xaxis.set_minor_locator(MultipleLocator(0.5))
        # Changing X-axis to match other figures
        ax.spines['bottom'].set_color('#B3D3BE')
        ax.spines['bottom'].set_linewidth(5)
        ax.tick_params(axis='x', which='both',
                       color='white', width=1,
                       direction='inout', length=5,
                       )
        ax.spines['bottom'].set_position(('outward', 5))
        ax.spines['bottom'].set_zorder(0)
        # Add ticks depending on graph position
        if i%ncols != 0:
            ax.tick_params(axis='y', which='both', left=False, labelleft=False)
        if i  < (df.shape[1] - ncols):
            ax.tick_params(axis='x', which='both', bottom=True, labelbottom=False)
    # Add figure-level labels
    fig.suptitle('Observed Recombination Rate on PAR1 per Donor')
    fig.supylabel('Rate (cM/Mb)')
    fig.supxlabel('Genomic Coordinates (Mb)')
    # Cleaning
    fig.tight_layout()
    # Plotting
    plt.savefig(path+'by_donor_PAR_recomb.svg',
                transparent=False,)
    plt.savefig(path+'by_donor_PAR_recomb.png',
            transparent=False,)
    return

def plot_averaged_physical_vs_genetic_distance(df, path, titletype):
    '''
    Plots averaged cumululative genetic distance across PAR1
    Inputs:
    df := Cumulative Recombination Frequency AUC dataframe
    path := directory folder to store images
    titletype := string to append to start (e.g. {proximal, distal})
    '''
    # Get limits normalized in megabases (Mb)
    start, end = 1/10**6, df.shape[0]/10**6
    fig, ax = plt.subplots(figsize=(10, 5), dpi=300)
    # Plotting average
    ax.plot((df.index/10**6)[telomeric_boundary:proposed_PAB], df.mean(axis=1)[telomeric_boundary:proposed_PAB],
            '-', c='#B3D4BE', label='Average (PAR1)',
            alpha=1.0,
            )
    ax.plot((df.index/10**6)[proposed_PAB:], df.mean(axis=1)[proposed_PAB:],
            '-', c='#74B088', label='Average (Contested PAB Interval)',
            alpha=1.0,
            )
    # Plotting proposed PAR boundary (PAB)
    plt.axvline(x=proposed_PAB/10**6, ymax=0.95,
                c='black', linestyle=(0, (1, 3, 6, 3, 1, 3)),
                linewidth=1.5,
                label='Proposed PAB',
                )
    # Plotting established PAR boundary (PAB)
    plt.axvline(x=established_PAB/10**6, ymax=0.95,
                c='black', linestyle=(0, (9, 3)),
                linewidth=1.5,
                label='Established PAB',
                )
    # Set graph limits
    ax.set_xlim(left=start, right=end+0.002)
    ax.set_ylim(bottom=0, top=50)
    # Annotate PAB
    ax.annotate('Proposed\nPAB', ha='center',
                xy=(proposed_PAB/10**6, 40), xytext=(proposed_PAB/10**6, 48),
                )
    ax.annotate('Established\nPAB', ha='center',
                xy=(established_PAB/10**6, 40), xytext=(established_PAB/10**6, 48),
                )
    # Labels
    plt.title('Averaged Genetic Distance on PAR1')
    ax.set_ylabel('Genetic Distance (cM)')
    ax.set_xlabel('Genomic Coordinates (Mb)')
    # Clean up graph
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_major_locator(MultipleLocator(10))
    # Changing X-axis to match other figures
    ax.spines['bottom'].set_color('#B3D3BE')
    ax.spines['bottom'].set_linewidth(5)
    ax.tick_params(axis='x', color='white', width=1,
                    direction='inout', length=5,
                    )
    ax.spines['bottom'].set_position(('outward', 5))
    ax.spines['bottom'].set_zorder(0)
    # Plotting
    fig.tight_layout()
    plt.savefig(path+'average_cumulative_PAR_recomb_'+titletype+'.svg',
                transparent=False,)
    plt.savefig(path+'average_cumulative_PAR_recomb_'+titletype+'.png',
            transparent=False,)
    return

def plot_donor_physical_vs_genetic_distance(df, path, titletype, ncols=2, figsize=(8, 10)):
    '''
    Plots genetic distance across PAR1 on a per donor basis
    Inputs:
    df := Cumulative Recombination Frequency AUC dataframe
    path := directory folder to store images
    titletype := string to append to start (e.g. {proximal, distal})
    ncols := number of columns to use in figure
    figsize := figure size in inches (tuple: width, height)
    '''
    # Get limits noramlized in megabases (Mb)
    start, end = 1/10**6, df.shape[0]/10**6
    # Create figure
    fig = plt.figure(figsize=figsize, dpi=300)
    # Define number of rows / columns for subplot grid
    nrows = df.shape[1] // ncols + (df.shape[1] % ncols > 0)
    for i, donor in enumerate(df.columns):
        ax = plt.subplot(nrows, ncols, (i + 1)) # new graph
        # Plotting
        ax.plot((df.index/10**6)[telomeric_boundary:proposed_PAB], df[donor][telomeric_boundary:proposed_PAB],
                '-', c='#B3D4BE', label=donor,
                alpha=1.0,
                )
        ax.plot((df.index/10**6)[proposed_PAB:], df[donor][proposed_PAB:],
                '-', c='#74B088', label=donor,
                alpha=1.0,
                )

        # Plotting proposed PAR boundary (PAB)
        plt.axvline(x=proposed_PAB/10**6, ymax=1.00,
                    c='black', linestyle=(0, (1, 3, 6, 3, 1, 3)),
                    linewidth=1.5,
                    label='Proposed PAB',
                    )
        # Plotting established PAR boundary (PAB)
        plt.axvline(x=established_PAB/10**6, ymax=1.00,
                    c='black', linestyle=(0, (9, 3)),
                    linewidth=1.5,
                    label='Established PAB',
                    )
        # Set graph limits
        ax.set_xlim(left=start, right=end+0.02)
        ax.set_ylim(bottom=0, top=50)
        # Labels
        ax.set_title((' ' + donor), loc='left', y=0.90, fontsize=12)
        # Clean up graph
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        # Set major tick locator for y-axis
        ax.yaxis.set_major_locator(MultipleLocator(5))
        ax.yaxis.set_major_formatter(FuncFormatter(format_ticks_for_donor_plots))
        # Set tick locators for x-axis
        ax.xaxis.set_major_locator(MultipleLocator(1.0))
        ax.xaxis.set_minor_locator(MultipleLocator(0.5))
        # Changing X-axis to match other figures
        ax.spines['bottom'].set_color('#B3D3BE')
        ax.spines['bottom'].set_linewidth(5)
        ax.tick_params(axis='x', which = 'both',
                       color='white', width=1,
                       direction='inout', length=5,
                       )
        ax.spines['bottom'].set_position(('outward', 5))
        ax.spines['bottom'].set_zorder(0)
        # Add ticks depending on graph position
        if i%ncols != 0:
            ax.tick_params(axis='y', which='both', left=False, labelleft=False)
        if i  < (df.shape[1] - ncols):
            ax.tick_params(axis='x', which='both', bottom=True, labelbottom=False)
    # Add figure-level labels
    fig.suptitle('Genetic Distance per Donor')
    fig.supylabel('Genetic Distance (cM)')
    fig.supxlabel('Genomic Coordinates (Mb)')
    # Cleaning
    fig.tight_layout()
    # Plotting
    plt.savefig(path+'by_donor_cumulative_PAR_recomb_'+titletype+'.svg',
                transparent=False,)
    plt.savefig(path+'by_donor_cumulative_PAR_recomb_'+titletype+'.png',
            transparent=False,)
    return

def AUC_histogram(df, PAB_AUC_df, path):
    '''
    Plots AUC across the entire PAR1 region
    or just the contested PAR1 boundary region.
    Inputs:
    df := pandas DataFrame, AUC across entire PAR1
    PAB_AUC_df := pandas DataFrame, AUC only in contested region on PAR1 boundary
    path := save directory
    '''
    # Plot AUC histogram
    max_AUC_array = [x for x in df.max()] # Getting maximums
    mean_AUC = np.mean(max_AUC_array)
    stddev_AUC = np.std(max_AUC_array)
    print('The average AUC and SD across the entire PAR1 is:')
    print('Mean AUC: {}\tSD: {}'.format(round(mean_AUC, 2), round(stddev_AUC, 2)))
    fig, ax = plt.subplots(figsize=(10,5), dpi=300)
    ax.hist(max_AUC_array, histtype='step',
            color='black',
            )
    # Add mean line
    plt.axvline(x=np.mean(max_AUC_array), ymin=0.00, ymax=0.15,
                c='black', linewidth=3,
                label='Mean')
    plt.axvline(x=np.mean(max_AUC_array), ymin=0.30, ymax=0.75,
                c='black', linewidth=3,
                label='Mean')
    # Annotate the mean
    x_limits = ax.get_xlim()
    x_range = x_limits[1] - x_limits[0]
    ax.annotate('Mean\nGenetic Distance', xy=((np.mean(max_AUC_array)-x_limits[0])/x_range, 0.20),
                xycoords='axes fraction', ha='center',
                )
    # Annotate with mean and standard deviation text
    plt.text(plt.xlim()[0] + 0.05 * (plt.xlim()[1] - plt.xlim()[0]),
             plt.ylim()[1] - 0.05 * (plt.ylim()[1] - plt.ylim()[0]),
             f'Mean ± Std Dev: {mean_AUC:.2f} ± {stddev_AUC:.2f}',
             va='top', ha='left', fontsize=10,
             bbox=dict(facecolor='white', alpha=0.5, boxstyle='round,pad=0.3'))
    # Clean up graph
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_title('Total Genetic Distance across the PAR1')
    ax.yaxis.set_major_locator(MultipleLocator(1))
    ax.set(xlabel='Total Genetic Distance (cM)', ylabel='Binned Donor Count')
    plt.tight_layout()
    # Plotting
    plt.savefig(path+'PAR_AUC_histogram.svg',
                transparent=False,)
    plt.savefig(path+'PAR_AUC_histogram.png',
            transparent=False,)

    # Plot AUC histogram for PAB interval
    max_PAB_AUC_array = [x for x in PAB_AUC_df.max()] # Getting maximums
    mean_AUC = np.mean(max_PAB_AUC_array)
    stddev_AUC = np.std(max_PAB_AUC_array)
    print('The average AUC and SD across the interval between the Proposed and Established PAB is:')
    print('Mean AUC: {}\tSD: {}'.format(round(mean_AUC, 2), round(stddev_AUC, 2)))
    fig, ax = plt.subplots(figsize=(10,5), dpi=300)
    ax.hist(max_PAB_AUC_array, histtype='step',
            color='black',
            )
    # Add mean line
    plt.axvline(x=np.mean(max_PAB_AUC_array), ymin=0.00, ymax=0.15,
                c='black', linewidth=3,
                label='Mean')
    plt.axvline(x=np.mean(max_PAB_AUC_array), ymin=0.30, ymax=0.75,
                c='black', linewidth=3,
                label='Mean')
    # Annotate the mean
    x_limits = ax.get_xlim()
    x_range = x_limits[1] - x_limits[0]
    ax.annotate('Mean\nGenetic Distance', xy=((np.mean(max_PAB_AUC_array)-x_limits[0])/x_range, 0.20),
                xycoords='axes fraction', ha='center',
                )
    # Annotate with mean and standard deviation text
    plt.text(plt.xlim()[0] + 0.05 * (plt.xlim()[1] - plt.xlim()[0]),
             plt.ylim()[1] - 0.05 * (plt.ylim()[1] - plt.ylim()[0]),
             f'Mean ± Std Dev: {mean_AUC:.2f} ± {stddev_AUC:.2f}',
             va='top', ha='left', fontsize=10,
             bbox=dict(facecolor='white', alpha=0.5, boxstyle='round,pad=0.3'))
    # Clean up graph
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_title('Total Genetic Distance in interval between Proposed and Established PAB')
    ax.yaxis.set_major_locator(MultipleLocator(1))
    ax.set(xlabel='Total Genetic Distance (cM)', ylabel='Binned Donor Count')
    plt.tight_layout()
    # Plotting
    plt.savefig(path+'PAB_AUC_histogram.svg',
                transparent=False,)
    plt.savefig(path+'PAB_AUC_histogram.png',
            transparent=False,)
    return

# Additional Functions

def genetic_distance_interval(df, start, stop, path, histogram_title, save_title):
    '''
    Given a set of genomic coordinates, calculate
    the genetic distance between them and return
    a histogram with the mean genetic distance.
    Inputs:
    df := Pandas DataFrame; recombinatino frequencies at each bp
    start := int; genomic coordinate marking start of interval
    stop := int; genomic coordinate marking end of interval
    '''
    # First, calculate exact genetic distance
    dx = 1/10**6 # Get interval
    # Define holder array for AUC
    AUC_array_cumulative = []
    # Cycle through donors
    for i, donor in enumerate(df.columns):
        # Calculating AUC only for bases in interval
        genomic_coordinates = (df.index/10**6)[start:stop]
        recomb_freq_interval = df[donor][start:stop]
        mask = ~np.isnan(recomb_freq_interval)
        genetic_distance = np.full(len(genomic_coordinates),
                                   fill_value=np.nan)
        genetic_distance[0] = 0 # Initialize count to 0
        # Iterate through recombination frequencies and calculate running AUC
        for h, j in enumerate(range(np.min(recomb_freq_interval.index), np.max(recomb_freq_interval.index)+1)):
            if mask[j] == True and h == 0:
                genetic_distance[h] = recomb_freq_interval[j] * dx
            elif mask[j] == True and h != 0:
                # Add the area of the rectangle formed by current RF
                # and the previous RF
                genetic_distance[h] = genetic_distance[h - 1] + (recomb_freq_interval[j] * dx)
            else:
                genetic_distance[h] = genetic_distance[h - 1]

        AUC_array_cumulative.append(genetic_distance)
    # Convert list of arrays to pandas DataFrame
    AUC_array_cumulative = pd.DataFrame(list(map(np.ravel, AUC_array_cumulative)))
    AUC_array_cumulative = AUC_array_cumulative.T
    AUC_array_cumulative.columns = df.columns

    # Plot histogram
    max_AUC_array = [x for x in AUC_array_cumulative.max()] # Getting maximums
    mean_AUC = np.mean(max_AUC_array)
    stddev_AUC = np.std(max_AUC_array)
    print('The average AUC and SD across the given interval is:')
    print('Mean AUC: {}\tSD: {}'.format(round(mean_AUC, 2), round(stddev_AUC, 2)))
    fig, ax = plt.subplots(figsize=(10,5), dpi=300)
    ax.hist(max_AUC_array, histtype='step',
            color='black',
            )
    # Add mean line
    plt.axvline(x=np.mean(max_AUC_array), ymin=0.00, ymax=0.15,
                c='black', linewidth=3,
                label='Mean')
    plt.axvline(x=np.mean(max_AUC_array), ymin=0.30, ymax=0.75,
                c='black', linewidth=3,
                label='Mean')
    # Annotate the mean
    x_limits = ax.get_xlim()
    x_range = x_limits[1] - x_limits[0]
    ax.annotate('Mean\nGenetic Distance', xy=((np.mean(max_AUC_array)-x_limits[0])/x_range, 0.20),
                xycoords='axes fraction', ha='center',
                )
    # Annotate with mean and standard deviation text
    plt.text(plt.xlim()[0] + 0.05 * (plt.xlim()[1] - plt.xlim()[0]),
             plt.ylim()[1] - 0.05 * (plt.ylim()[1] - plt.ylim()[0]),
             f'Mean ± Std Dev: {mean_AUC:.2f} ± {stddev_AUC:.2f}',
             va='top', ha='left', fontsize=10,
             bbox=dict(facecolor='white', alpha=0.5, boxstyle='round,pad=0.3'))
    # Clean up graph
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_title(histogram_title)
    ax.yaxis.set_major_locator(MultipleLocator(1))
    ax.set(xlabel='Total Genetic Distance (cM)', ylabel='Binned Donor Count')
    plt.tight_layout()
    # Plotting
    plt.savefig(path+'AUC_histogram_'+save_title+'.svg',
                transparent=False,)
    plt.savefig(path+'AUC_histogram_'+save_title+'.png',
            transparent=False,)
    return

if __name__ == '__main__':
    # Get input arguments
    parser = argparse.ArgumentParser(description='Find recombination files and plot results.')
    parser.add_argument('-d', '--dir', help='Path to parent directory to grab files. Defaults to current folder.',
                        default='./')
    parser.add_argument('-o', '--out', help='Path to directory to save. Defaults to current folder.',
                        default='./')
    parser.add_argument('-c', '--cols', help='Number of columns for donor separated graphs. Defaults to 2.',
                        default=2, type=int)
    parser.add_argument('-w', '--width', help='Figure width of donor separated graphs. Defaults (8) inches.',
                        default=8, type=int)
    parser.add_argument('-l', '--length', help='Figure length (height) of donor separated graphs. Defaults (10) inches.',
                        default=10, type=int)
    parser.add_argument('-f', '--font', help='Path to .ttf or .ttc font file (set up for Helvetica). Defaults to DejaVu Sans.',
                        default='', type=str)
    args = vars(parser.parse_args())

    dir_path = args['dir']
    print('File Directory:', dir_path)
    save_dir = args['out']
    print('Save Directory:', save_dir)
    cols = args['cols']
    figsize = (args['width'], args['length'])
    print('Donor Separated Graphs has Figure Size: {} (inches) and {} columns.'.format(figsize, cols))
    
    font_path = args['font']
    if font_path != '':
        fm.fontManager.addfont(font_path)
        plt.rcParams['font.family'] = 'sans-serif'
        plt.rcParams['font.sans-serif'] = ['Helvetica']
        # Set the svg.fonttype parameter to 'none' to ensure text remains as text in the SVG
        plt.rcParams['svg.fonttype'] = 'none'
    else:
        plt.rcParams['svg.fonttype'] = 'none'

   # Import all relevant txt files to build pandas DataFrame
    print('Finding files with recombination rates...')
    files = find_files('NC*.bedGraph*', dir_path)
    print('Found {} files.'.format(len(files)))

    # For each file, build DataFrame
    print('Building Recombination Frequency DataFrame...')
    # Get recombination frequencies (RF) in list of numpy arrays
    holder_array = []
    global_min, global_max = 1, 1 # 1-indexed
    for file in files:
        df = pd.read_csv(file, sep='\t',
                         skiprows=1, header=None,
                         )
        df = df.rename(columns={0:'Chromosome', 1:'Start', 2:'End', 3:'RF'})
        # Get minimum and global maximum
        if df['Start'].min() <= global_min:
            global_min = df['Start'].min()
        if df['End'].max() >= global_max:
            global_max = df['End'].max()
        # Append df to list
        holder_array.append(get_recomb_array(df))
    # Get donor names
    donor_names = []
    for file in files:
        donor_names.append(re.findall(r'NC\d+', file)[0])
    # Get number of sperm used
    num_sperm = []
    for file in files:
        num_sperm.append(int(re.search(r'across (\d+) Sperm', pd.read_csv(file, sep='\t', nrows=0).columns[0]).group(1)))
    total_sperm = np.sum(num_sperm)
    # Combine together
    donor_names = ['{} (n={})'.format(nc, n) for nc, n in zip(donor_names, num_sperm)]
    # Removing "raggedness" from arrays
    df = np.full((len(holder_array), global_max-global_min), fill_value=np.NaN)
    for donor, array in enumerate(df):
        for i in range(len(holder_array[donor])):
            array[i] = holder_array[donor][i]
    # Making pandas DataFrame
    df = pd.DataFrame.from_dict(dict(zip(donor_names, df)))
    # Reordering DataFrame lexicographically
    df = df[sorted(donor_names, key=custom_key)]
    print('Built Recombination Frequency Dataframe.')

    # Getting recombination frequency AUC dataframes
    print('Calculating proximal Recombination Frequency AUC wrt telomere...')
    AUC_df_proximal, PAB_AUC_df_proximal = get_cum_AUC_proximal_from_telomere(df=df)
    print('Calculating distal Recombination Frequency AUC wrt est. PAR1 boundary...')
    AUC_df_distal, PAB_AUC_df_distal = get_cum_AUC_distal_from_PAR(df=df)

    # Sanity Check readouts
    print('Checking that proximal and distal calculations are (near) equivalent...')
    check_AUC(AUC_df_proximal, AUC_df_distal, portion_of_par='entirety')
    check_AUC(PAB_AUC_df_proximal, PAB_AUC_df_distal, portion_of_par='contested boundary region')

    # Plot recombination frequency at each point
    print('Plotting average recombination frequency in PAR1...')
    plot_averaged_rf(df, path=save_dir)
    print('Plotting donor-level recombination frequency in PAR1...')
    plot_donor_rf(df, path=save_dir, ncols=cols, figsize=figsize)

    # Plot cumulative recombination frequency
    print('Plotting average cumulative recombination frequency in PAR1...')
    plot_averaged_physical_vs_genetic_distance(df=AUC_df_proximal, path=save_dir, titletype='proximal')
    plot_averaged_physical_vs_genetic_distance(df=AUC_df_distal, path=save_dir, titletype='distal')
    print('Plotting donor-level cumulative recombination frequency in PAR1...')
    plot_donor_physical_vs_genetic_distance(df=AUC_df_proximal, path=save_dir, titletype='proximal',
                                            ncols=cols, figsize=figsize)
    plot_donor_physical_vs_genetic_distance(df=AUC_df_distal, path=save_dir, titletype='distal',
                                            ncols=cols, figsize=figsize)

    # Plot histogram of total AUCs
    print('Plotting histogram of total AUCs...')
    AUC_histogram(df=AUC_df_distal, PAB_AUC_df=PAB_AUC_df_distal, path=save_dir)

    # Calculate genetic differences for any other intervals of interest
    # and return histogram plot of maximum genetic distance per donor
    print('Calculating genetic distance for Hinch et al marker interval...')
    genetic_distance_interval(df=df,
                              start=hinch_proximal_marker,
                              stop=hinch_distal_marker,
                              path=save_dir,
                              histogram_title='Genetic Distance in interval between markers used in Hinch et al.',
                              save_title='Hinch_markers',)

    print('Done.')
