#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 14:56:05 2020

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

cann_pop = sim.Population(n, sim.IF_cond_alpha(**cell_params),
                         label = "cann_pop")

inhib_pop = sim.Population(1, sim.IF_cond_alpha(**cell_params),
                         label = "inhib_pop")

spike_source = sim.Population(1, sim.SpikeSourceArray(
    spike_times = spike_times))

spike_2_conn = sim.Projection(spike_source, cann_pop[spiky:spiky+1], sim.AllToAllConnector(),
                      sim.StaticSynapse(weight = 0.002, delay = 0.1))

cann_2_cann = sim.Projection(cann_pop, cann_pop, sim.FromListConnector(cann_connector, column_names=["weight"]),
                          sim.StaticSynapse(weight = 0.0001, delay = 75))

cann_2_inh = sim.Projection(cann_pop, inhib_pop, sim.AllToAllConnector(),
                            sim.StaticSynapse(weight = 0.02, delay = 0.1),
                            receptor_type = "excitatory")

inh_2_cann = sim.Projection(inhib_pop, cann_pop, sim.AllToAllConnector(),
                            sim.StaticSynapse(weight = 0.2, delay = 0.1),
                            receptor_type = "inhibitory")

spike_source.record('spikes')
cann_pop.record(('v','spikes'))
inhib_pop.record(('v','spikes'))

sim.run(5000.0)

from pyNN.utility.plotting import Figure, Panel
cann_data = cann_pop.get_data().segments[0]
cann_vm = cann_data.filter(name="v")[0]
inh_data = inhib_pop.get_data().segments[0]
inh_vm = inh_data.filter(name="v")[0]

Figure(
       Panel(cann_vm, ylabel = "Cann Membrane potential (mV)", yticks = True),
       Panel(cann_data.spiketrains, yticks = True),
       Panel(inh_vm, ylabel = "Inh Membrane potential (mV)", yticks = True),
       Panel(inh_data.spiketrains, xlabel = "Time (ms)", yticks = True, xticks = True)
)

sim.end()
