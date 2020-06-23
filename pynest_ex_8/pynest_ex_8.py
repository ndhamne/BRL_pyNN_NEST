#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 07:00:33 2020

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

spike_time_for_id_1 = [100.]
spike_time_for_id_2 = [500., 505., 510., 515.]
spike_time_for_id_3 = [1200., 1205., 1210., 1215.]
spike_time_for_id_4 = [2000., 2005., 2010., 2015.]
spike_time_for_id_9 = [3000., 3005., 3010., 3015.]
#[3000., 3003., 3005., 3008., 3010., 3013., 3015., 3017., 3020.]

cann_pop = sim.Population(n, sim.IF_cond_alpha(**cell_params),
                         label = "cann_pop")

inhib_pop = sim.Population(1, sim.IF_cond_alpha(**cell_params),
                         label = "inhib_pop")

spike_source_1 = sim.Population(1, sim.SpikeSourceArray(
    spike_times = spike_time_for_id_1))
spike_source_2 = sim.Population(1, sim.SpikeSourceArray(
    spike_times = spike_time_for_id_2))
spike_source_3 = sim.Population(1, sim.SpikeSourceArray(
    spike_times = spike_time_for_id_3))
spike_source_4 = sim.Population(1, sim.SpikeSourceArray(
    spike_times = spike_time_for_id_4))
spike_source_9 = sim.Population(1, sim.SpikeSourceArray(
    spike_times = spike_time_for_id_9))

spike_2_conn_for_1 = sim.Projection(spike_source_1, cann_pop[1:2], sim.OneToOneConnector(),
                      sim.StaticSynapse(weight = 0.002, delay = 0.1))
spike_2_conn_for_2 = sim.Projection(spike_source_2, cann_pop[2:3], sim.OneToOneConnector(),
                      sim.StaticSynapse(weight = 0.002, delay = 0.1))
spike_2_conn_for_3 = sim.Projection(spike_source_3, cann_pop[3:4], sim.OneToOneConnector(),
                      sim.StaticSynapse(weight = 0.002, delay = 0.1))
spike_2_conn_for_4 = sim.Projection(spike_source_4, cann_pop[4:5], sim.OneToOneConnector(),
                      sim.StaticSynapse(weight = 0.002, delay = 0.1))
spike_2_conn_for_9 = sim.Projection(spike_source_9, cann_pop[9:10], sim.OneToOneConnector(),
                      sim.StaticSynapse(weight = 0.002, delay = 0.1))



cann_2_cann = sim.Projection(cann_pop, cann_pop, sim.FromListConnector(cann_connector, column_names=["weight"]),
                          sim.StaticSynapse(weight = 0.0001, delay = 75))

cann_2_inh = sim.Projection(cann_pop, inhib_pop, sim.AllToAllConnector(),
                            sim.StaticSynapse(weight = 0.02, delay = 0.1),
                            receptor_type = "excitatory")

inh_2_cann = sim.Projection(inhib_pop, cann_pop, sim.AllToAllConnector(),
                            sim.StaticSynapse(weight = 0.2, delay = 0.1),
                            receptor_type = "inhibitory")

#spike_source.record('spikes')
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