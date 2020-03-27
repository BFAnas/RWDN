"""
    This script generates random water distribution networks as INP files.
"""

from RWDN_core import RandomWaterDistributionNetwork
import pandas as pd
import osmnx as ox
import networkx as nx
import random
import sys
import os.path

# gmap_key = input("Enter gmap_key (Google maps API): \n")
gmap_key = "AIzaSyC83ArsL5MYzYP91KV13WvAP2Dmix8M-xU"

# Creating the ouput directory
output_dir = input("Enter output path: \n")
if os.path.exists(output_dir) == False:
    output_dir = os.getcwd() + '/WDN_output'
try:
    os.mkdir(output_dir)
except OSError:
    print ("Creation of the directory %s already exists" % output_dir)
else:
    print ("Successfully created the directory %s " % output_dir)


try:
    number_of_WDN = int(input("Enter the number of Water Distribution Networks to be generated: \n"))
except ValueError:
    print("This is not a whole number.")


# Loading cities data from which the location of the streets will be generated
world_cities = pd.read_csv("/home/anas/Desktop/Thesis/simplemaps_worldcities_basicv1.6/worldcities.csv")

# the minimal population for the cities to be used for the generation of the Virtual Water Network
min_population = 100000

# Average demand values in gal/day/acre from Water Distribution Systems Handbook - McGraw Hill (2000)
demand_values = {'Low_density_residential': 1670, 'Med_density_residential': 2610, 'High_density_residential': 4160, 'Single_family_residential ': 2300, 'Multifamily residential': 4160, 'Office_commercial': 2030, 'Retail_commercial': 2040, 'Light_industrial': 1620, 'Heavy_industrial': 2270, 'Parks': 2020, 'Schools': 1700}

# Average hourly demand values in m3/h/hectare
demand_values = {key: v*0.00039 for key,v in demand_values.items()}

# Roughness
roughness_values = [100, 110, 130]

# Reservoir's elevation over the highest point
reservoir_heads = [70, 80]

# Possible pipes diameters values in m
# pipes_diameters = [0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.13, 0.15, 0.18, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]

pipes_diameters = [0.08, 0.1, 0.13, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]


for i in range(number_of_WDN):

    while True:
        try:
            # distance = random.randint(500, 5000)
            distance = 1000

            rd = RandomWaterDistributionNetwork(gmap_key, world_cities, min_population, distance, demand_values, roughness_values, reservoir_heads, pipes_diameters)

            # Generate the layout from street network data 
            while True:
                try:
                    G = rd.generate_random_graph()
                except G.size() < 200:
                    continue
                break

            print ('Layout n: {} ' .format(i) + 'was generated successfully')
            # Plot the generated network
            # fig, ax = ox.plot_graph(G, fig_height=8, node_size=20)

            # Add elevation values to the nodes
            rd.add_elevation(G)

            # Add node demands 
            rd.add_node_demands(G)

            # Convert G to undirected graph
            G = nx.Graph(G)

            # Add reservoirs as nodes and connecting them to the graph
            rd.add_reservoirs(G)

            # Add edge roughness
            rd. add_edge_roughness(G)

            # Create WaterNetworkModel object for wntr
            wn = rd.create_wn(G)

            # Run simulation for pipe sizing
            wn = rd.pipe_sizing(wn)

            # Divide the network to sectors and add valves between sectors
            wn = rd.add_valves(wn)

        except rd.stats(wn)['Mean_pressure'] < 20 or rd.stats(wn)['Mean_pressure'] > 100:
            continue
        break

    # Save the WDN as inp file
    filename = os.path.join(output_dir, 'WDN{}' .format(i) + '.inp')
    rd.save_wn(wn, filename)

print("Successfully created {} " .format(number_of_WDN) + " virtual WDNs located in " + output_dir)
