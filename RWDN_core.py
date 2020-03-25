import os
import numpy as np
import math
import osmnx as ox
import networkx as nx
import pandas as pd
import random
import matplotlib.pyplot as plt
import community as com
from networkx.algorithms import community
import wntr


class RandomWaterDistributionNetwork:

    def __init__(self, gmap_key, world_cities, min_population, distance, demand_values, roughness_values, reservoir_heads, pipes_diameters):

        '''
        This class is for constructing Random Virtual Water Distribution Networks from OpenStreetMap streets layout using osmnx module, and using Wntr module for diameter sizing   


        gmap_key(str): key for google maps API for downloading the elevation values for the nodes

        world_cities(pd.DataFrame): pandas DataFrame from the free csv file in: https://simplemaps.com/data/world-cities

        min_population(int): the minimal population for the cities to be used for the generation of the Virtual Water Network

        distance(int): is the distance of the nodes from the center point that are going to be kept for constructing the graph

        demand_values(dict): a dictionnary where the keys are the type of soil occupation, and values are peak hour demands in m3/h/hectar

        roughness_values(list): is a list of roughness values that will be assigned randomly to pipes

        reservoir_heads(list): is a list of heads that will be assigned randomly to the reservoirs 

        number_of_reservoirs(int)

        pipe_diameters(list): list of possible diameters of pipes in the networks (in m)
        '''

        self.gmap_key = gmap_key
        self.world_cities = world_cities
        self.min_population = min_population
        self.distance = distance
        self.demand_values = demand_values
        self.roughness_values = roughness_values
        self.reservoir_heads = reservoir_heads
        self.pipes_diameters = pipes_diameters

        # A list to strore problematic pipes (pipes where the start node is also the end node)
        self.problematic_pipes = []
        self.number_of_reservoirs = 2


    def generate_random_graph(self):
        '''
        Returns a random graph using the osmnx module
        '''
        big_cities = self.world_cities[self.world_cities.population > self.min_population].reset_index().drop(columns=['index'])
        random_index = random.randint(0, len(big_cities))

        random_city_coord = (big_cities.loc[random_index, 'lat'], big_cities.loc[random_index, 'lng'])

        r1 = random.randint(-100, 100)/5000
        r2 = random.randint(-100, 100)/5000
        random_point = (random_city_coord[0]+r1, random_city_coord[1]+r2)

        return ox.graph_from_point(random_point, distance=self.distance, network_type='drive')

    def add_elevation(self, G):
        '''
        Add elevation attribute to the nodes of the input graph
        '''
        G = ox.add_node_elevations(G, api_key=self.gmap_key)


    def add_node_demands(self, G, area = 1):
        '''
        Add demand attribute to the nodes of the input graph. The graph is divided to communities and each community is assigned a demand value randomly from demand_values

        area(float): The average area in hectare covered by each of the nodes
        '''
        G_copy = nx.Graph(G)
        n_nodes = len(G)
        n_communities = round(n_nodes/4)
        community_generator = community.asyn_fluidc(G_copy, n_communities)
        demands = {}

        # Converting demand values to peak demands in m3/s
        demand_values = {key: 4*area*v/3600 for key,v in self.demand_values.items()}           

        for i in range(n_communities):
            # Construct a dictionary of demands for each community
            # (community nodes are the keys and demand values are the values)
            demands_list = list(demand_values.values())
            demands.update(dict.fromkeys(next(community_generator), {'demand' : demands_list[i%len(demand_values)]}))

        # Add the demand as attributes to the nodes in the graph
        nx.set_node_attributes(G, demands) 

    def add_edge_roughness(self, G):
        '''
        Adds the roughness as attribute to the edges in the graph
        '''
        nx.set_edge_attributes(G, random.choice(self.roughness_values), 'roughness') 

    def add_reservoirs(self, G):
        '''
        Adds number_of_reservoirs reservoirs to the graph G
        '''
        G_copy = nx.Graph(G)
        n_edges = G_copy.number_of_edges()
        self.number_of_reservoirs = min(round(n_edges/1000) +1, round(self.distance/500) +1)
        n_communities = self.number_of_reservoirs
        community_generator = community.asyn_fluidc(G_copy, n_communities)
        reservoir_nodes =[]
        for i in range(n_communities):
            # Make a subgraph from the community nodes
            subG = G.subgraph(next(community_generator))
            # Get the maximum elevation node from the subgraph
            max_elevation_node = max(dict(subG.nodes).items(), key=lambda x: x[1]['elevation'])[0]
            Reservoir_attr = {'x': G.nodes[max_elevation_node]['x']+ random.randint(-9, 9)/10000, 'y': G.nodes[max_elevation_node]['y'] + random.randint(-9, 9)/10000, 'elevation': random.choice(self.reservoir_heads) + G.nodes[max_elevation_node]['elevation'], 'demand': 0}
            G.add_node('Reservoir{}' .format(i), **Reservoir_attr)
            G.add_edge('Reservoir{}' .format(i), max_elevation_node)
            G['Reservoir{}'.format(i)][max_elevation_node]['length'] = 50
            G['Reservoir{}'.format(i)][max_elevation_node]['roughness'] = random.choice(self.roughness_values)

    def create_wn(self, G):
        '''
        Creates wntr.network.WaterNetworkModel from the graph G
        '''
        wn = wntr.network.WaterNetworkModel()
        G_junctions = G.copy()
        for i in range(self.number_of_reservoirs):
            G_junctions.remove_node('Reservoir{}' .format(i))

        n_nodes = len(G_junctions)
        n_edges = G.number_of_edges()

        nodes_list = list(G_junctions.nodes(data=True))
        edges_list = list(G.edges(data=True))

        wn.add_pattern('pat1', [0.56, 0.37, 0.28, 0.33, 0.45, 0.67, 1.13, 2.02, 1.78, 1.6, 1.56, 1.45, 1.18, 1.19, 1.07, 0.97, 0.67, 0.9, 0.82, 1.22, 1.57, 1.3, 1.05, 0.92])

        for i in range(n_nodes):
            node = nodes_list[i]
            wn.add_junction('{}'.format(node[0]), base_demand=node[1]['demand'], elevation=node[1]['elevation'], coordinates=(node[1]['x'], node[1]['y']))


        for i in range(self.number_of_reservoirs):
            reservoir = G.nodes['Reservoir{}' .format(i)]
            wn.add_reservoir('Reservoir{}' .format(i), base_head=reservoir['elevation'], coordinates=(reservoir['x'], reservoir['y']))
        

        for i in range(n_edges):
            edge = edges_list[i]
            if edge[0] == edge[1] or edge[2]['length'] > 2000: self.problematic_pipes.append(i)
            else:
                wn.add_pipe('{}' .format(i), edge[0], edge[1], length=edge[2]['length'], diameter=0.2, roughness=edge[2]['roughness'],
                minor_loss=0.0, status='OPEN')

        return wn


    def pipe_sizing(self, wn):
        '''
        Adjusting the diameters of the pipes to get realistic values for the pressures and velocities in the network
        '''
        # wn = wntr.network.WaterNetworkModel(inp_file)
        # Simulate hydraulics
        sim = wntr.sim.EpanetSimulator(wn)
        results = sim.run_sim()
        velocity = results.link['velocity'].values[0]
        pipes_list = set(range(len(velocity)))-set(self.problematic_pipes)

        def bigger_diameter(diameter):
            i = 0
            while diameter >= self.pipes_diameters[i] and i < len(self.pipes_diameters)-1 : i+=1
            return self.pipes_diameters[i]

        def smaller_diameter(diameter):
            i = len(self.pipes_diameters)-1
            while diameter <= self.pipes_diameters[i] and i > 0: i-=1
            return self.pipes_diameters[i]

        sim = wntr.sim.EpanetSimulator(wn)
        results = sim.run_sim()
        pressure = results.node['pressure'].values[0][:-2].mean()
        flow = results.link['flowrate'].values[0]
        # diam = abs(4*flow/3.14)
        c = 0
        counter = len(pipes_list)
        starting_diam = self.pipes_diameters[0]
        while pressure < 0 :
            starting_diam = bigger_diameter(starting_diam)
            for link in pipes_list:
                pipe = wn.get_link('{}'.format(link))
                pipe.diameter = bigger_diameter(starting_diam)
            sim = wntr.sim.EpanetSimulator(wn)
            results = sim.run_sim()
            pressure = results.node['pressure'].values[0][:-2].mean()
        # while counter != 0 and c < len(self.pipes_diameters):
        while  c < len(self.pipes_diameters) and counter > len(pipes_list)/200:
            counter = 0
            c += 1
            sim = wntr.sim.EpanetSimulator(wn)
            results = sim.run_sim()
            velocity = results.link['velocity'].values[0]
            print(c)
            for link in pipes_list:
                pipe = wn.get_link('{}'.format(link))
                node1 = pipe.start_node_name
                node2 = pipe.end_node_name
                pressure1 = results.node['pressure'][node1][0]
                pressure2 = results.node['pressure'][node2][0]
                pressure = (pressure1 + pressure2)/2
                print(pressure)
                if velocity[link] < 0.5 and pressure > 70:
                    pipe.diameter = smaller_diameter(pipe.diameter)
                    print(pipe.diameter, c)
                    counter += 1
                if velocity[link] > 1.5 and pressure < 40:
                    pipe.diameter = bigger_diameter(pipe.diameter)
                    print(pipe.diameter, c)
                    counter += 1
                if velocity[link] > 2.5:
                    pipe.diameter = bigger_diameter(pipe.diameter)
                    print(pipe.diameter, c)
                    counter += 1
                if velocity[link] < 0.01:
                    pipe.diameter = smaller_diameter(pipe.diameter)
                    print(pipe.diameter, c)
                    counter += 1
        return wn

    def add_valves(self, wn, plot_sect='False'):
        '''
        Divide the network to sectors using Louvain-Algorithm, and add valves between sectors. If plot_sect = 'True' it also plots the network with different colors for the sectors
        '''
        G = wn.get_graph()
        G = nx.Graph(G)
        # community.generate_dendrogram(G)
        sectors = com.best_partition(G)
        counter = 0
        pipes_list = set(range(wn.num_links))-set(self.problematic_pipes)
        for link in pipes_list:
            pipe = wn.get_link('{}'.format(link))
            node1 = pipe.start_node_name
            node2 = pipe.end_node_name
            sect1 = sectors['{}'.format(node1)]
            sect2 = sectors['{}'.format(node2)]
            if sect1 != sect2 :
                counter += 1
                diameter = pipe.diameter
                print(diameter)
                wn.add_valve('S{}'.format(counter), node1, node2, diameter, 'TCV', 0)
        
        if plot_sect == 'True':
            nodes, edges = wntr.graphics.plot_network(wn, node_attribute=sectors)

        return wn

            
    def save_wn(self, wn, filename):
        ''' 
        Save the network as inp file
        '''
        # Change the base_demand back to average demand instead of peak demand, and introduce demand paterns
        for node_name in wn.node_name_list[:-self.number_of_reservoirs]:
            node = wn.get_node(node_name)
            node.demand_timeseries_list[0].base_value = 0.25 * node.base_demand
            node.pattern = 'pat1'

        wn.write_inpfile(filename, units='CMH')

    @property
    def stats(self, wn):
        '''
        Calculates some stats of the network
        '''

        sim = wntr.sim.EpanetSimulator(wn)
        results = sim.run_sim()

        # Mean velocity:
        mean_velocity = results.link['velocity'].values[0].mean()

        # Min velocity:
        min_velocity = results.link['velocity'].values[0].min()

        # Max_velocity:
        max_velocity = results.link['velocity'].values[0].max()

        # Mean pressure:
        mean_pressure = results.node['pressure'].values[0][:-2].mean()

        # Min pressure: 
        min_pressure = results.node['pressure'].values[0][:-2].min()

        # Max pressure:
        max_pressure = results.node['pressure'].values[0][:-2].max()

        #

        return {'Mean velocity' : mean_velocity, 'Min velocity' : min_velocity, 'Max_velocity' : max_velocity, 'Mean pressure' : mean_pressure, 'Min pressure' : min_pressure, 'Max pressure' : max_pressure}






