#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 12:25:36 2020

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

cann_pop = sim.Population(n, sim.IF_cond_alpha(**cell_params),
                         label = "cann_pop")
dist = np.zeros((n, n))
weights = np.zeros((n, n))
cann_connector = list()

for i in range(n):
    for j in range(n):
        dist[i][j] = abs(i - j)%(n)
        if dist[i][j] > n/2:
            dist[i][j] = n - dist[i][j]
        weights[i][j] = round(weight_to_spike * \
                    (1 / (sig * np.sqrt(2 * np.pi)) * \
                    (np.exp(-np.power(dist[i][j] - mu, 2.) \
                    / (2 * np.power(sig, 2.))))), 2)
        cann_connector.append((i, j, weights[i][j]))#, delay_cann2cann))
print("Weight matrix:\n", weights)                    

spike_times = [1000., 2000.]

spike_source = sim.Population(1, sim.SpikeSourceArray(
    spike_times = spike_times))

conn = sim.Projection(spike_source, cann_pop[spiky:spiky+1], sim.AllToAllConnector(),
                      sim.StaticSynapse(weight = 0.002, delay = 0.1))

n2n_conn = sim.Projection(cann_pop, cann_pop, sim.FromListConnector(cann_connector, column_names=["weight"]),
                          sim.StaticSynapse(weight = 0.0001, delay = 100))

spike_source.record('spikes')
cann_pop.record(('v','spikes'))

sim.run(5000.0)

from pyNN.utility.plotting import Figure, Panel
data1 = cann_pop.get_data().segments[0]
vm = data1.filter(name="v")[0]

Figure(
       Panel(vm, ylabel = "Membrane potential (mV)", yticks = True),
       Panel(data1.spiketrains, xlabel = "Time (ms)", yticks = True, xticks = True)
)

sim.end()
