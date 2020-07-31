#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 29 17:49:11 2020

@author: nikhildhamne
"""

import pyNN.nest as sim
import numpy as np

n = 10
spiky = n//2

mu = 0.
sig = 1

weight_to_spike = 0.025 #0.01 #0.05
delay_cann2cann = 1.

cell_params = {
    'tau_m'      : 20.0,   # (ms)
    'tau_syn_E'  : 2.0,    # (ms)
    'tau_syn_I'  : 4.0,    # (ms)
    'e_rev_E'    : 0.0,    # (mV)
    'e_rev_I'    : -70.0,  # (mV)
    'tau_refrac' : 2.0,    # (ms)
    'v_rest'     : -60.0,  # (mV)
    'v_reset'    : -70.0,  # (mV) 
    'v_thresh'   : -60.0,  # (mV)
    'cm'         : 0.5}
sim.setup()

dist = np.zeros((n, n))
weights = np.zeros((n, n))
hd_cann_connector = list()
pos_rot_connector = list()
neg_rot_connector = list()
hd_conj_connector = list()

for i in range(n):
    for j in range(n):
        dist[i][j] = abs(i - j)%(n)
        if dist[i][j] > n/2:
            dist[i][j] = n - dist[i][j]
        weights[i][j] = round(weight_to_spike * \
                    (1 / (sig * np.sqrt(2 * np.pi)) * \
                    (np.exp(-np.power(dist[i][j] - mu, 2.) \
                    / (2 * np.power(sig, 2.))))), 2)
        hd_cann_connector.append((i, j, weights[i][j]))#, delay_cann2cann))
#print("Weight matrix:\n", weights)                    


hd_spike_times = [2000., 2005., 2010., 2015.]

spikes_in = 20 #5 #10 #50
pos_spike_times = [None] * spikes_in
for i,j in zip(range(spikes_in), range(4000, 4100, 2)):
    pos_spike_times[i] = j
print(pos_spike_times)
#pos_spike_times = [2000., 2005., 2010., 2015.]
neg_spike_times = []#[4000., 4005., 4010., 4015.]

# Populations

hd_cann_pop = sim.Population(n, sim.IF_cond_alpha(**cell_params),
                         label = "hd_cann_pop")

inhib_pop = sim.Population(n, sim.IF_cond_alpha(**cell_params),
                         label = "inhib_pop")

#pos_rot_conj = sim.Population(n, sim.IF_cond_alpha(**cell_params),
pos_rot_conj = sim.Population(n, sim.IF_cond_alpha(),
                          label = "pos_rot_conj")

#neg_rot_conj = sim.Population(n, sim.IF_cond_alpha(**cell_params),
neg_rot_conj = sim.Population(n, sim.IF_cond_alpha(),
                          label = "neg_rot_conj")

# Spike Sources

hd_spike_src = sim.Population(1, sim.SpikeSourceArray(
    spike_times = hd_spike_times))

pos_spike_src = sim.Population(1, sim.SpikeSourceArray(
    spike_times = pos_spike_times))

neg_spike_src = sim.Population(1, sim.SpikeSourceArray(
    spike_times = neg_spike_times))

# Connections between spike sources and ppopulation 

spike_2_hd = sim.Projection(hd_spike_src, hd_cann_pop[3:4], sim.OneToOneConnector(),
                      sim.StaticSynapse(weight = 0.002, delay = 0.1))

pos_spikes = sim.Projection(pos_spike_src, pos_rot_conj, sim.AllToAllConnector(),
                            sim.StaticSynapse(weight = 0.2, delay = 0.1))

neg_spikes = sim.Projection(neg_spike_src, neg_rot_conj, sim.AllToAllConnector(),
                            sim.StaticSynapse(weight = 0.2, delay = 0.1))

# Connections/ Projections between main populations

hd_cann_2_hd_cann = sim.Projection(hd_cann_pop, hd_cann_pop, 
                                   sim.FromListConnector(hd_cann_connector, column_names=["weight"]),
                                   sim.StaticSynapse(weight = 0.0001, delay = 75))

hd_cann_2_inh = sim.Projection(hd_cann_pop, inhib_pop, sim.AllToAllConnector(),
                               sim.StaticSynapse(weight = 0.02, delay = 0.1),
                               receptor_type = "excitatory")

inh_2_hd_cann = sim.Projection(inhib_pop, hd_cann_pop, sim.AllToAllConnector(),
                               sim.StaticSynapse(weight = 0.5, delay = 0.1),
                               receptor_type = "inhibitory")

for i in range(n):
    pos_rot_connector.append((i, (i+1)%10, 0.002))
#print(pos_rot_connector)   
pos_rot_conj_2_hd_cann = sim.Projection(pos_rot_conj, hd_cann_pop, 
                                        sim.FromListConnector(pos_rot_connector, column_names=["weight"]),
                                        receptor_type = "excitatory")

for i in range(n):
    neg_rot_connector.append((i, (i-1)%10, 0.002))
    
neg_rot_conj_2_hd_cann = sim.Projection(neg_rot_conj, hd_cann_pop, 
                                        sim.FromListConnector(neg_rot_connector, column_names=["weight"]),
                                        sim.StaticSynapse(delay = 100),
                                        receptor_type = "excitatory")
for i in range(n):
    hd_conj_connector.append((i, i, weight_to_spike))
    
hd_cann_2_pos_rot = sim.Projection(hd_cann_pop, pos_rot_conj,
                                   sim.FromListConnector(hd_conj_connector, column_names=["weight"]),
                                   sim.StaticSynapse(delay = 100),
                                   receptor_type = "excitatory")

hd_cann_2_neg_rot = sim.Projection(hd_cann_pop, neg_rot_conj,
                                   sim.FromListConnector(hd_conj_connector, column_names=["weight"]),
                                   receptor_type = "excitatory")

#spike_source.record('spikes')
hd_cann_pop.record(('v','spikes'))
inhib_pop.record(('v','spikes'))
pos_rot_conj.record(('v','spikes'))
neg_rot_conj.record(('v','spikes'))

sim.run(10000.0)

from pyNN.utility.plotting import Figure, Panel
cann_data = hd_cann_pop.get_data().segments[0]
cann_vm = cann_data.filter(name="v")[0]
inh_data = inhib_pop.get_data().segments[0]
inh_vm = inh_data.filter(name="v")[0]

pos_rot_data = pos_rot_conj.get_data().segments[0]
pos_rot_vm = pos_rot_data.filter(name="v")[0]
neg_rot_data = neg_rot_conj.get_data().segments[0]
neg_rot_vm = neg_rot_data.filter(name="v")[0]


Figure(
       #Panel(cann_vm, ylabel = "Cann Membrane potential (mV)", yticks = True),
       Panel(cann_data.spiketrains, ylabel = "HD spikes", yticks = True),
       #Panel(inh_vm, ylabel = "Inh Membrane potential (mV)", yticks = True),
       Panel(inh_data.spiketrains, ylabel = "Inhib. spikes", yticks = True),
       Panel(pos_rot_data.spiketrains, ylabel = "Pos Rot spikes", yticks = True),
       Panel(neg_rot_data.spiketrains, xlabel = "Time (ms)", ylabel = "Neg Rot spikes", yticks = True, xticks = True)
       #Panel(pos_rot_vm, ylabel = "Pos Rot pot (mV)", yticks = True),
       #Panel(neg_rot_vm, xlabel = "Time (ms)", ylabel = "Neg Rot pot (mV)", yticks = True, xticks = True)
)

sim.end()