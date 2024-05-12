#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 11 19:13:14 2024

@author: dilettaolliaro
"""

from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib
import csv
import re


# Figure Settings

fsize = 150
legend_size = 200
label_size = 220
title_size = 195
tuplesize = (90, 70)
line_size = 25
tick_size = 180
l_pad = 40


# Parameters
N = 200
T = 100

mus = 1
mub = 0.25

ps = 0.8
pb = 0.2

filename = f'Results/experiment_overLoad_N{N}_T{T}_ps{ps:.2f}_mu_s{mus:.2f}_mu_b{mub:.2f}.csv'  
lab = fr'$p_s$ = {ps}, $\mu_s$ = {mus:.2f}, $\mu_b$ = {mub:.2f}'

load = []    
respTime = []    

with open(filename) as csv_file:
    
    df = pd.read_csv(filename, delimiter = ',')

    for index, row in df.iterrows():
         
         load.append(float(row['Load']))
         respTime.append(float(row['RespTime Total']))
             

############################## OVERALL RESP TIME ##############################      
    
plt.figure(dpi=1200)
plt.rc('font',**{'family':'serif','serif':['Palatino']})
plt.rc('text', usetex=True)
matplotlib.rcParams['font.size'] = fsize
matplotlib.rcParams['xtick.major.pad'] = 8
matplotlib.rcParams['ytick.major.pad'] = 8
fix, ax = plt.subplots(figsize=tuplesize)          

ax.plot(load, respTime, color = 'crimson', label = lab, ls = 'solid', lw = line_size)
plt.axvline(x = 0.76, color = 'black', label = 'stability', ls = 'dotted', lw = line_size)

ax.set_xlabel("Load", fontsize=label_size)
ax.set_ylabel("Avg. Overall Resp. Time $\quad[$s$]$", fontsize=label_size)
ax.set_title(f"Avg. Overall Resp. Time vs. Load, N = {N}, T = {T}", fontsize=title_size)
plt.xscale('linear')
plt.yscale('linear')
ax.tick_params(axis='both', which='major', labelsize=tick_size, pad = l_pad)
ax.tick_params(axis='both', which='minor', labelsize=tick_size, pad = l_pad)
ax.legend(fontsize = legend_size)
ax.grid()
plt.savefig(f'Plots/respTimeVsLoad_N{N}_T{T}_ps{ps:.2f}_mu_s{mus:.2f}_mu_b{mub:.2f}.pdf')
