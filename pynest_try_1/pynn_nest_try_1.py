#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 27 15:37:00 2020

@author: nikhildhamne
"""

import pyNN.nest as sim
from pyNN.utility.plotting import Figure, Panel

# variables names and values
n_neurons = 10
sim_time = 6000.0

# Build the network
sim.setup()

# spike source
#generate_spike_times = [[0], [500], [1000], [1500], [2000], [2500], [3000], 
#                       [3500], [4000], [4500]] #, [5000]]
#spike_source = sim.Population(n_neurons, 
#                              sim.SpikeSourceArray(
#                                  spike_times = generate_spike_times))

spike_source = sim.Population(n_neurons, 
                              sim.SpikeSourceArray(
                                  spike_times = [0, 1000, 3000]))


# populations and projections
layer_1 = sim.Population(n_neurons, sim.IF_curr_alpha(),
                         label = "layer_1")

layer_2 = sim.Population(n_neurons, sim.IF_curr_alpha(),
                         label = "layer_2")

proj_src_2_l1 = sim.Projection(spike_source, layer_1,
                               sim.OneToOneConnector())

proj_l1_2_l2 = sim.Projection(layer_1, layer_2, 
                              sim.OneToOneConnector())

sim.run(sim_time)

layer_1.record(('spikes', 'v'))
layer_2.record(['spikes', 'v'])
#layer_1.record(['v', 'spikes'])
#layer_2.record(['v', 'spikes'])

# plot the spikes
data_1 = layer_1.get_data().segments[0]
vm = data_1.filter(name = "v")
Figure(
       Panel(vm, ylabel = "Memb. Pot."),
       Panel(data_1.spiketrains, xlabel = "Time (ms)", xticks = True)
)

sim.end()
