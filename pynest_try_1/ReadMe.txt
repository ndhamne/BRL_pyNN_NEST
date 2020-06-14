Here, PyNN is used as front end to create the network and NEST is used as a backend.
Start by importing necessary dependencies, and declare variables. For simplicity, the network has only 10 neurons. 

A spike source - SpikeSourceArray to inject spikes into the network “layer1” at particular times listed in the array “spike_times”.

After creating the spike_source and the layer1 population of neurons, the spike source is connected to the layer1. The inbuilt “OneToOneconnector” is used to connect every one neuron in spike source to layer 1 populations. 

We then start recording the membrane potential and spikes as the simulation beings. Then they are plotted against time.
