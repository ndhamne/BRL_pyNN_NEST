#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 12:54:50 2020

@author: nikhildhamne
"""
import pyNN.nest as sim

n = 5
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

neuron1 = sim.Population(n, sim.IF_cond_alpha(**cell_params),
                         label = "neuron1")

neuron2 = sim.Population(n, sim.IF_cond_alpha(**cell_params),
                         label = "neuron2")

generate_spike_times = [0., 1020., 1040., 1060., 
                        1080., 1100., 1120., 1140.,
                        1160., 2000.]

spike_source = sim.Population(n, sim.SpikeSourceArray(
    spike_times = generate_spike_times))

conn = sim.Projection(spike_source, neuron1, sim.AllToAllConnector(),
                      sim.StaticSynapse(weight = 0.002, delay = 1.))

n2n_conn = sim.Projection(neuron1, neuron2, sim.OneToOneConnector(),
                          sim.StaticSynapse(weight = 0.002, delay = 1.))

spike_source.record('spikes')
neuron1.record(('v','spikes'))
neuron2.record(('v','spikes'))

sim.run(5000.0)

#print neuron1.get_spike_counts()
from pyNN.utility.plotting import Figure, Panel
data1 = neuron1.get_data().segments[0]
data2 = neuron1.get_data().segments[0]
vm = data1.filter(name="v")[0]

Figure(
       Panel(vm, ylabel = "Membrane potential (mV)", yticks = True),
       Panel(data1.spiketrains, xlabel = "Time (ms)", xticks = True),
       Panel(data2.spiketrains, xlabel = "Time (ms)", xticks = True)
)

sim.end()
