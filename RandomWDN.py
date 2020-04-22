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
import community as com
import matplotlib as plt
import math
import wntr
import csv

# gmap_key = input("Enter gmap_key (Google maps API): \n")
gmap_key = "Insert your Gmap key"

# Creating the ouput directory
output_dir = os.getcwd() + '/WDN_output'

number_of_WDN = 100

# Loading cities data from which the location of the streets will be generated
world_cities = pd.read_csv(os.getcwd() + '/worldcities.csv')

# the minimal population for the cities to be used for the generation of the Virtual Water Network
min_population = 100000

# Average demand values in gal/day/acre from Water Distribution Systems Handbook - McGraw Hill (2000)
demand_values = {'Low_density_residential': 1670, 'Med_density_residential': 2610, 'High_density_residential': 4160, 'Single_family_residential ': 2300, 'Multifamily residential': 4160, 'Office_commercial': 2030, 'Retail_commercial': 2040, 'Light_industrial': 1620, 'Heavy_industrial': 2270, 'Parks': 2020, 'Schools': 1700}

# Average hourly demand values in m3/h/hectare
demand_values = {key: v*0.00039 for key,v in demand_values.items()}

# Roughness
roughness_values = [100, 110, 130]

# Reservoir's elevation over the highest point
reservoir_heads = [30, 40, 50]

# Possible pipes diameters values in m
pipes_diameters = [0.04, 0.06, 0.08, 0.1, 0.13, 0.15, 0.18, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]

# pipes_diameters = [0.08, 0.1, 0.11, 0.13, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]


j=0
while j < number_of_WDN:
# try:
    # distance = random.randint(1000, 3000)
    distance = 2000

    rd = RandomWaterDistributionNetwork(gmap_key, world_cities, min_population, distance, demand_values, roughness_values, reservoir_heads, pipes_diameters)

    # Generate the layout from street network data 
    while True:
        try:
            G = rd.generate_random_graph()
            rd.add_elevation(G)
            if len(G) < 200: 
                print ('Number of nodes < 200, this graph will be discarded') 
                continue
            else: break
        except Exception:
            continue
    print('Graph n: {} generated' .format(j))
    # Add elevation values to the nodes 
    G = ox.project_graph(G)
    print('Graph n: {} projected UTM-WGS84' .format(j))
    G = nx.Graph(G)

    # Clean the graph from clustered nodes, self loops and parallel edges
    G = rd.clean_graph(G, 15)
    print('Graph n: {} cleaned from too closely clustered nodes' .format(j))

    # Add node demands 
    rd.add_node_demands(G)
    print('Demand values added to nodes' .format(j))

    # Create main distribution network
    # subG = rd.generate_main_distr(G)
    subG = rd.main_network(G)
    subG = rd.clean_cycles(subG)
    print('Main distribution network created')

    # Add reservoirs as nodes and connecting them to the graph
    subG = rd.add_reservoirs(subG)
    # print(len(subG))
    print('Reservoirs added')
    F = nx.compose(G, subG)
    # F = rd.clean_graph(F, 0)
    # Add edge roughness
    rd.add_edge_roughness(F)
    print('Roughness added to the Graph {}'.format(j))

    # Create WaterNetworkModel object for wntr
    wn = rd.create_wn(F, subG)

    # Run simulation for pipe sizing
    wn = rd.pipe_sizing(wn)
    print('Pipe sizing performed')

    # Divide the network to sectors and add valves between sectors
    wn = rd.add_valves(wn, subG)
    stats = rd.stats(wn)
    print('Valves have been added to the network')
    
    if min(stats['Pressure']) > 0:
        filename = os.path.join(output_dir, 'WDN{}' .format(j) + '.inp')
        rd.save_wn(wn, filename)
        # wn = wntr.network.WaterNetworkModel(filename)
        print('File n: {} ' .format(j) + 'was generated successfully')

        print('Main valves list:')
        print(rd.main_valves)
        
        out = os.path.join(output_dir, 'stats/WDN_out{}' .format(j) + '.csv')
        w = csv.writer(open(out, "w"))
        for key, val in stats.items():
            w.writerow([key, val])
        
        j+=1

    else:
        continue
    # except Exception:
        # continue

