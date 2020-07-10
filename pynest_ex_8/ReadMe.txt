In this simulation the adjacent neurons are spiked at different times during the simulation to see if the bump moves smoothly which seems to be the case as seen from generated plot.

A neuron which is not in the adjacent neighbourhood is spiked with similar number of spikes to see if the activity shifts abruptly, but it doesn't move. This is neuro-plausible as well becasue the agent would only move with smoothly as opposed to the kidnapped robot scenario.
NOTE : The commented spike times for a farther neuron. If between to non-adjacent neurons, a neuron is injected with spikes sufficiently based on the confidence level of the location, the activity does shift in such a scenario with continues activity of the stable bump as well. The situation has two competing hypothesis foa a location. However, this requires very strong stimulation of that farther neuron. The activity of that case is also plotted and uploded to this commit. 

To get 8a_addl.png, comment the working spike times and uncomment the ones below it.
