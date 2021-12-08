import collections
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import numpy as np


def get_summary(graph):
    degree_sum = sum([d for n, d in graph.degree()])
    degree_centrality_sum = sum(nx.degree_centrality(graph).values())
    closeness_centrality_sum = sum(nx.closeness_centrality(graph).values())
    betweenness_centrality_sum = sum(nx.betweenness_centrality(graph).values())
    number_of_nodes = graph.number_of_nodes()
    shannon_entropy = calculate_shannon_entropy(graph)
    wiener_index = nx.wiener_index(graph)

    d = {
        'Number of Nodes': [number_of_nodes],
        'Number of Edges': [graph.number_of_edges()],
        'Diameter': [nx.diameter(graph) if nx.algorithms.components.is_connected(graph) else 'graph is not connected'],
        'Average Degree': degree_sum / number_of_nodes,
        'Average Degree Centrality': degree_centrality_sum / number_of_nodes,
        'Average Closeness Centrality': closeness_centrality_sum / number_of_nodes,
        'Average Betweenness Centrality': betweenness_centrality_sum / number_of_nodes,
        'Average Clustering Coefficient': [nx.average_clustering(graph)],
        'Shannon Entropy': shannon_entropy,
        'Wiener Index:': f'{wiener_index:.5f}'
    }

    df = pd.DataFrame(data=d)
    df = df.T
    print(df)


def get_giant_components(graph, number):
    gcc = sorted(nx.connected_components(graph), key=len, reverse=True)
    graphs = []
    for x in range(number):
        g0 = graph.subgraph(gcc[x])
        graphs.append(g0)

    return graphs


def calculate_shannon_entropy(graph):
    degrees = graph.degree()
    df = pd.DataFrame(degrees)
    v = len(df.index)
    grouped_degrees = df.groupby([1]).size()

    h = 0
    for di in grouped_degrees:
        pi = di / v
        h += pi * np.log10(pi)

    h *= -1

    return h


def plot_degree_histogram(graph):
    degree_sequence = sorted([d for n, d in graph.degree()], reverse=True)  # degree sequence
    degree_count = collections.Counter(degree_sequence)
    deg, cnt = zip(*degree_count.items())
    plot_histogram(deg, cnt, "Degree Histogram", "Count", "Degree")


def plot_histogram(values, cnt, tile, xlabel, ylabel):
    fig, ax = plt.subplots()
    plt.bar(values, cnt, width=0.80, color='b')

    plt.title(tile)
    plt.ylabel(xlabel)
    plt.xlabel(ylabel)
    ax.set_xticks([d + 0.4 for d in values])
    ax.set_xticklabels(values)

    plt.show()
