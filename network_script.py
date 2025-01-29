# %%
from pyneuroml import pynml
import urllib.request, json 
#import requests
import os
from neuroml import *
from neuroml.utils import component_factory, validate_neuroml2
from pyneuroml import pynml
from pyneuroml.lems import LEMSSimulation
import neuroml.writers as writers
import random
from pyneuroml.pynml import read_lems_file, read_neuroml2_file, write_lems_file, write_neuroml2_file
from pathlib import Path

# %%
'''https://pyneuroml.readthedocs.io/en/development/pyneuroml.io.html
pyneuroml.io.read_lems_file
pyneuroml.io.read_neuroml2_file
pyneuroml.io.write_lems_file
pyneuroml.io.write_neuroml2_file'''

# use the above functions to read LEMS and NeuroML files
print(Path(os.getcwd()).parent)
cell_data_dir = os.path.join(Path(os.getcwd()).parent, "NMC_Model/NMC.NeuronML2/")
network_data_dir = os.path.join(Path(os.getcwd()).parent, "NMC_Model/")  # contains JSONs for network and connectivity data

files = os.listdir(cell_data_dir)
cell_files = [f for f in files if f.endswith("cell.nml")]  # .cell.nml -- contains arborization location details
net_files = [f for f in files if f.endswith("net.nml")]  # .net.nml -- contains location (0,0,0) in x-y-z space
current_files = [f for f in files if (f.endswith(".nml") and 'current' in f)]  # .current_clamp.nml -- contains current clamp details
print(cell_files[:5])
print(net_files[:5])
print(current_files[:5])

# validate the NeuroML files
'''for file in cell_files:
    print("Validating: " + file)
    validate_neuroml2(cell_data_dir + file)
for file in net_files:
    validate_neuroml2(cell_data_dir + file)
for file in current_files:
    validate_neuroml2(cell_data_dir + file)
'''

# %%
# read in the network and connection properties from jsons: circuit.json, layers.json
with open(network_data_dir + 'circuit.json') as f:
    circuit_data = json.load(f)  # layer thickness and neuron counts
with open(network_data_dir + 'layers.json') as f:
    layers_data = json.load(f)  # layer-specific me-type counts
with open(network_data_dir + 'pathways_anatomy.json') as f:
    pathways_anatomy_data = json.load(f)  # contains the synapse counts/distributions for each m_type->m_type connection
with open(network_data_dir + 'pathways_physiology.json') as f:
    pathways_physiology_data = json.load(f)  # contains the synapse properties for each m_type->m_type connection
layers_data

# %%
nml_doc = component_factory("NeuroMLDocument", id="Cortex_Network")

### Create the network
net = nml_doc.add("Network", id="Cortex_Network", validate=False)
net.type="networkWithTemperature" 
net.temperature="34.0degC"

neuron_population_dict = {}  # neuron_population_dict[m_type][e_type] = population_object

for layer in layers_data:
    layer_dict = layers_data[layer]
    e_type_counts = layer_dict['No. of neurons per electrical types']
    m_type_counts = layer_dict['No. of neurons per morphological types']

    for m_type in m_type_counts:
        m_type_cells = [f for f in cell_files if m_type in f]
        m_count = m_type_counts[m_type]

        # for the e-types that match this m-type, build a relative probability distribution of these e-types based on their counts
        me_type_probs = {}
        e_count_sum = 0

        for e_type in e_type_counts:
            me_type_matches = [f for f in m_type_cells if e_type in f]
            for match in me_type_matches:
                me_type_probs[match] = e_type_counts[e_type] / len(me_type_matches)
            if len(me_type_matches) > 0:
                e_count_sum += e_type_counts[e_type]
        
        # normalize the probabilities
        for e_type in me_type_probs:
            me_type_probs[e_type] /= e_count_sum
        #print(me_type_probs)

        # probabilities should sum to 1
        if abs(1 - sum(me_type_probs.values())) > 0.0001 and len(me_type_probs.keys()) > 0:
            print("Error: probabilities do not sum to 1")
            print(sum(me_type_probs.values()))
            print(me_type_probs)

        # sample from the me-type distribution to get the number of cells of each type
        me_type_counts = {}  # maps me-type to # neurons of that type
        for i_morp in range(m_count):
            m_type_cells = list(me_type_probs.keys())
            m_type_probs = list(me_type_probs.values())
            if len(m_type_cells) == 0 or len(m_type_probs) == 0:
                continue
            me_type_sample = random.choices(m_type_cells, m_type_probs)[0]
            if me_type_sample not in me_type_counts:
                me_type_counts[me_type_sample] = 1
            else:
                me_type_counts[me_type_sample] += 1

        # create the neurons and place them in the network for this layer. Record in neuron_population_dict
        for me_type_sample in me_type_counts:
            print("Creating " + str(me_type_counts[me_type_sample]) + " cells of type " + me_type_sample)
            # read in the cell file
            cell_model = read_neuroml2_file(cell_data_dir + me_type_sample)

            # create the neuron population object and place them in the network for this layer
            size_pop1 = me_type_counts[me_type_sample]
            nml_doc.add("IncludeType", href=cell_data_dir + me_type_sample)

            population_obj = component_factory("Population", id="Exc", component=cell_model.id, size=size_pop1, type="population")
            # Set optional color property. Note: used later when generating plots
            ##pop0.add("Property", tag="color", value="0 0 .8")
            net.add(population_obj, force=True)

            # add to neuron_population_dict
            if m_type not in neuron_population_dict:
                neuron_population_dict[m_type] = {}
            neuron_population_dict[m_type][me_type_sample] = population_obj


# %%
total_connections = 0
for layer in layers_data:
    layer_dict = layers_data[layer]
    e_type_counts = layer_dict['No. of neurons per electrical types']
    m_type_counts = layer_dict['No. of neurons per morphological types']
    # now go through the connections and create the synapses
    for i_m1 in range(len(m_type_counts)):
        m1 = list(m_type_counts.keys())[i_m1]
        for i_m2 in range(len(m_type_counts)):
            m2 = list(m_type_counts.keys())[i_m2]
            m1m2_key = m1 + ":" + m2
            if m1m2_key in pathways_anatomy_data:
                print("Creating synapses " + m1m2_key) 
                syn_anatomy = pathways_anatomy_data[m1m2_key]
                syn_physiology = pathways_physiology_data[m1m2_key]

                # count the number of me-type pairs
                num_me_type_pairs = m_type_counts[m1] * m_type_counts[m2]
                
                # number of synapses per me-type pair
                num_me_synapses = syn_anatomy["total_synapse_count"] / num_me_type_pairs
                connection_prob = syn_anatomy["connection_probability"]
                mean_num_syn_per_connection = syn_anatomy["mean_number_of_synapse_per_connection"]
                std_num_syn_per_connection = syn_anatomy["number_of_synapse_per_connection_std"]

                # load physiology data
                g_syn = syn_physiology["gsyn_mean"]
                decay_mean = syn_physiology["decay_mean"]
                synapse_type = syn_physiology["synapse_type"]
                ei_type, stp_type = synapse_type.split(", ")


                # simple mono-exponential synapse
                # <Property name="weight" dimension="none" defaultValue="1"/>

                # <Parameter name="tauDecay" dimension="time" description="Time course of decay"/>
                # <Parameter name="gbase" dimension="conductance" description="Baseline conductance, generally the maximum conductance following a single spike"/>
                # <Parameter name="erev" dimension="voltage" description="Reversal potential of the synapse"/>

                syn0 = nml_doc.add(
                    "ExpOneSynapse", id=m1m2_key.replace(":", "_"), gbase="65nS", erev="0mV", tau_decay= str(decay_mean) + "ms"
                )
                #print("me-types:")
                #print(neuron_population_dict[m1])
                #print(neuron_population_dict[m2])
                # for every pair of me-types, create a fraction of the synapses
                for me1 in neuron_population_dict[m1]:
                    me_pop1 = neuron_population_dict[m1][me1]
                    for me2 in neuron_population_dict[m2]:
                        me_pop2 = neuron_population_dict[m2][me2]

                        # Create connections and inputs
                        proj_count = 0

                        projection = Projection(
                                            id="Proj_",
                                            presynaptic_population=me_pop1.id,
                                            postsynaptic_population=me_pop2.id,
                                            synapse=syn0.id,
                                        )

                        net.projections.append(projection)
                        for i in range(me_pop1.size):
                            for j in range(me_pop2.size):
                                if random.random() <= connection_prob: # probablistic connection...
                                    connection = ConnectionWD(
                                        id=proj_count,
                                        pre_cell_id="%s[%i]" % (me_pop1.id, i),
                                        post_cell_id="%s[%i]" % (me_pop2.id, j),
                                        weight=random.random(),
                                        delay='0ms'
                                    )
                                    projection.add(connection)
                                    proj_count += 1
                        print("\tAdded %i connections from %s to %s" % (proj_count, me1, me2))
                        total_connections += proj_count
            print("Total so far: %i" % total_connections + "\n\n")



# %%
# stimulus: L4 spiny stellate neurons (simulating thalamic input) as a pulse generator at 20ms, duration 0.2ms, amplitude on 10^-1 nA scale
for me_pop in neuron_population_dict["L4_SS"]:
    for i in range(0, me_pop.size):
        # pulse generator as explicit stimulus
        pg = nml_doc.add(
            "PulseGenerator",
            id="pg_exc_%i" % i,
            delay="20ms",
            duration="0.2ms",
            amplitude="%f nA" % (0.3 + 0.1 * random.random()),
        )
        
        exp_input = net.add(
            "ExplicitInput", target="%s[%i]" % (me_pop.id, i), input=pg.id
        )

print(nml_doc.summary())

nml_net_file = 'Cortex_Network.net.nml'
writers.NeuroMLWriter.write(nml_doc, nml_net_file)

print("Written network file to: " + nml_net_file)
pynml.validate_neuroml2(nml_net_file)


