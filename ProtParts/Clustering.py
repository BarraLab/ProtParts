import networkx as nx
import operator
from sklearn.metrics import silhouette_samples
import numpy as np
from .utils import hobohm1

class Cluster:

    """
    Cluster object
    """

    def __init__(self, clusters={}):
        """
        Parameters
        ----------
        cluster : dict
            Clustered sequences
        """
        self.clusters = clusters
    

    def __len__(self):
        """
        Returns
        -------
        length : int
            Number of clusters
        """
        return len(self.clusters)
    

    def __eq__(self, other):
        """
        Parameters
        ----------
        other : Cluster
            Other clusters object
        
        Returns
        -------
        result : bool
            True if two clusters are identical, False otherwise
        """
        if len(self.clusters) != len(other.clusters):
            return False
        
        clusters_1_set = {frozenset(clusters) for clusters in self.clusters.values()}
        clusters_2_set = {frozenset(clusters) for clusters in other.clusters.values()}
        
        return clusters_1_set == clusters_2_set


    def __ne__(self, other):
        """
        Parameters
        ----------
        other : Cluster
            Other clusters object
        
        Returns
        -------
        result : bool
            True if two clusters are not identical, False otherwise
        """
        return not self.__eq__(other)
    

    def items(self):
        """
        Returns
        -------
        items : list
            List of clusters items
        """
        return self.clusters.items()


    def num_data(self, by='sum'):
        """
        Parameters
        ----------
        by: str/int
            'max': maximum number of data among clusters
            'min': minimum number of data among clusters
            'sum': sum of number of data in all clusters
            str/int: number of data in the cluster with the given id
        
        Returns
        -------
        size : int
            Number of sequences in all clusters
        """
        size_list = [len(c) for c in self.clusters.values()]

        if by == 'max':
            size = max(size_list)
        elif by == 'min':
            size = min(size_list)
        elif by == 'sum':
            size = sum(size_list)
        elif isinstance(by, int) or isinstance(by, str):
            size = len(self.clusters[by])
        
        return size
    

    def sort(self, reverse=False):
        """
        Sort clusters by size

        Parameters
        ----------
        reverse : bool
            Reverse sort
        
        Returns
        -------
        result : Cluster
            Sorted clusters
        """
        result = Cluster(dict(sorted(self.clusters.items(), key=lambda x:len(x[1]), reverse=reverse)))
        return result
 
    
    # def list_all(self):
    #     """
    #     Returns
    #     -------
    #     result : list
    #         List of all data
    #     """
    #     result = []
    #     for c in self.clusters.values():
    #         result.extend(c)
    #     return result
    

    def index(self, data=None):
        """
        Parameters
        ----------
        data : str
            Sequence id
        
        Returns
        -------
        idx : int
            Cluster index of the data
        """
        inverted = {v:k for k, values in self.clusters.items() for v in values}
        if data is None:
            idx = inverted
        elif data not in inverted:
            raise ValueError(f"Invalid data: {data}")
        else:
            idx = inverted[data]
        return idx



    # def evaluate_matrix(self, method=None, matrix=None):
    #     """
    #     Evaluate the cluster based on similarity/distance matrix

    #     Parameters
    #     ----------
    #     method : str
    #         Evaluation method
    #     matrix : np.array
    #         Similarity/Distance matrix
        
    #     Returns
    #     -------
    #     metric : float
    #         Evaluation result
    #     """

    #     if method == 'sihouette':
    #         return self._silhouette(matrix)

    
    def silhouette(self, measurement):
        """
        Silhouette score

        Parameters
        ----------
        measurement : tuple
            List of measurement (seq1, seq2, measurement)
        
        Returns
        -------
        metric : float
            Silhouette score
        """
        data_list = list(self.index().keys())
        data_label = list(self.index().values())

        if len(self.clusters) == 1 or len(self.clusters) == len(data_list):
            return None, (data_list, data_label, None)
        else:
            pivot = np.ones((len(data_list), len(data_list))) * 11
            data_dict = dict(zip(data_list, range(len(data_list))))
            for row in measurement:
                seq1, seq2, measure = row[0], row[1], row[2]
                if (seq1 in data_dict) and (seq2 in data_dict):
                    i = data_dict[seq1]
                    j = data_dict[seq2]
                    pivot[i, j] = measure
            
            np.fill_diagonal(pivot, 0)
            sample_silhouette_values = silhouette_samples(pivot, data_label, metric='precomputed')
            metric = np.mean(sample_silhouette_values)

            return metric, (data_list, data_label, sample_silhouette_values)


class Clustering:

    """
    Clustering methods
    """

    def __init__(self, threshold, method, measurement_type='distance'):
        """
        Parameters
        ----------
        threshold : float
            Threshold for clustering
        method : str
            Clustering method
        measurement_type : str
            Type of measurement (distance or similarity)
        
        Returns
        -------
        None

        Raises
        ------
        ValueError
            If measurement_type is not 'distance' or 'similarity'
        """
        self.threshold = threshold
        #self.num_clusters = num_clusters

        if method not in ['graph', 'hobohm1']:
            raise ValueError('Invalid clustering method: {}'.format(method))
        self.method = method

        if measurement_type not in ['distance', 'similarity']:
            raise ValueError('Invalid measurement type: {}'.format(measurement_type))
        self.measurement_type = measurement_type
    

    def clustering(self, sequences, measurement):
        """
        Sequence clustering based on the measurement

        Parameters
        ----------
        sequences : dict
            Dict of sequences
        measurement : tuple
            List of measurement (seq1, seq2, measurement)
        
        Returns
        -------
        result : Cluster
            Clustered sequences
        """
        if self.measurement_type == 'distance':
            op = operator.le
        elif self.measurement_type == 'similarity':
            op = operator.ge

        if self.method == 'graph':
            G = self._graph(sequences, measurement, op)
            result = {idx:sorted(list(component)) for idx, component in enumerate(nx.connected_components(G))}
            clusters = Cluster(result)
        elif self.method == 'hobohm1':
            clusters = self._hobohm1(sequences, measurement, op)
        
        return clusters


    def _graph(self, sequences, measurement, op):
        """
        Create a graph from the sequences and measurement

        Parameters
        ----------
        sequences : dict
            Dict of sequences
        measurement : tuple
            List of measurement (seq1, seq2, measurement)

        Returns
        -------
        result : Cluster
            Clustered sequences
        """
        G = nx.Graph()
        #nodes = list(map(lambda x:x.id, sequences))
        nodes = list(sequences.keys())
        G.add_nodes_from(nodes)

        # filter out self-self measurement and based on threshold
        #measurement = list(filter(lambda x:(x[0] != x[1]) and (x[0] in sequences) and (x[1] in sequences) and (op(x[2], self.threshold)), measurement))
        # select the first three columns of the measurement
        measurement = [(x[0], x[1], x[2]) for x in measurement if (x[0] != x[1]) and (x[0] in sequences) and (x[1] in sequences) and (op(x[2], self.threshold))]
        
        # add edges
        G.add_weighted_edges_from(measurement)

        return G
    

    def _hobohm1(self, sequences, measurement, op):
        """
        Hobohm1 clustering algorithm

        Parameters
        ----------
        sequences : dict
            Dict of sequences
        measurement : tuple
            List of measurement (seq1, seq2, measurement)
        
        Returns
        -------
        result : Cluster
            Clustered sequences
        """
        
        result = hobohm1(sequences, measurement, self.threshold, op, reduce_redundancy=False)

        return Cluster(result)


    def _avg_distance_to_neighbors(self, G, node):
        """
        Calculate the average distance of a node to its neighbors

        Parameters
        ----------
        node : str
            Node id

        Returns
        -------
        avg_dist : float
            Average distance to neighbors
        """

        neighbors = list(G.neighbors(node))
        if len(neighbors) == 0:
            return float('inf')  # Return infinity if a node has no neighbors
        distances = [G[node][nbr].get('weight', 1) for nbr in neighbors]  # Default to 1 if unweighted
        return sum(distances) / len(distances)


    def _ratio_centrality(self, G):
        """
        Calculate the ratio centrality of each node

        Returns
        -------
        ratios : dict
            Dict of ratio centrality of each node
        """

        ratios = {}
        for node in G.nodes():
            avg_dist = self._avg_distance_to_neighbors(G, node)
            # Calculate the average distance of each neighbor to their neighbors
            neighbor_dists = [self._avg_distance_to_neighbors(G, nbr) for nbr in G.neighbors(node)]
            # Avoid division by zero
            if neighbor_dists:
                avg_neighbor_dist = sum(neighbor_dists) / len(neighbor_dists)
                ratios[node] = avg_dist / avg_neighbor_dist if avg_neighbor_dist else float('inf')
            else:
                ratios[node] = float('inf')  # No neighbors case
        return ratios


    def _optimize_ratio(self, sequences, measurement):
        """
        Optimize the cluster based on ratio centrality

        Parameters
        ----------
        sequences : dict
            Dict of sequences
        measurement : tuple
            List of measurement (seq1, seq2, measurement)
        
        Returns
        -------
        result : Cluster
            Optimized cluster
        """
        
        #clust = Clustering(threshold=threshold, method='graph', measurement_type='distance')
        G = self._graph(sequences, measurement, operator.le)
        cluster = Cluster({idx:sorted(list(component)) for idx, component in enumerate(nx.connected_components(G))})
        silhouette_score, silhouette_score_samples = cluster.silhouette(measurement)

        silhouette_score_tmp = 99
        silhouette_score_current = silhouette_score
        best_step = {'graph':G, 'cluster':cluster, 'silhouette_score':silhouette_score, 'silhouette_score_samples':silhouette_score_samples}
        
        iter = 0
        #print(f"Iter\tCluster\tNode\tSilhouette_iter\tSilhouette_highest\tSilhouette_current")
        while silhouette_score_tmp > silhouette_score_current:
            iter += 1
            #print(f"Iter {iter}")
            #print([len(c) for c in best_step['silhouette_score_samples']])
            negative_clusters = {c for (i, c, v) in zip(*best_step['silhouette_score_samples']) if v < 0}
            #print(negative_clusters)
            
            best_step_tmp = {'silhouette_score':-1}
            silhouette_score_tmp = best_step['silhouette_score']
            silhouette_score_current = best_step['silhouette_score']

            for c in negative_clusters:
                elements = best_step['cluster'].clusters[c]
                G_tmp = best_step['graph'].subgraph(elements)

                ratios = self._ratio_centrality(G_tmp)
                nodes_tmp = sorted(ratios, key=ratios.get, reverse=True)
                
                """
                # remove the node with the highest ratio centrality
                for node in nodes_tmp:
                    sequence_new = {k:v for k, v in sequences.items() if (k != node) and (k in best_step['cluster'].index())}
                    measurement_new = list(filter(lambda x:(x[0] != node) and (x[1] != node), measurement))
                    G_new = self._graph(sequence_new, measurement_new, operator.le)
                    cluster_new = Cluster({idx:sorted(list(component)) for idx, component in enumerate(nx.connected_components(G_new))})
                    silhouette_score_new, silhouette_score_samples_new = cluster_new.silhouette(measurement_new)
                    
                    print(f"{iter}\t{c}\t{node}\t{silhouette_score_new:.5f}\t{silhouette_score_tmp:.5f}\t{best_step['silhouette_score']:.5f}\t{silhouette_score_new > silhouette_score_tmp}\t{best_step_tmp['silhouette_score'] > best_step['silhouette_score']}")
                    
                    # best_step_tmp record the best step for each negative cluster
                    if silhouette_score_new > silhouette_score_tmp:
                        best_step_tmp['graph'] = G_new
                        best_step_tmp['cluster'] = cluster_new
                        best_step_tmp['silhouette_score'] = silhouette_score_new
                        best_step_tmp['silhouette_score_samples'] = silhouette_score_samples_new
                        silhouette_score_tmp = silhouette_score_new
                    else:
                        break
                """
                node = nodes_tmp[0]
                sequence_new = {k:v for k, v in sequences.items() if (k != node) and (k in best_step['cluster'].index())}
                measurement_new = list(filter(lambda x:(x[0] != node) and (x[1] != node), measurement))
                G_new = self._graph(sequence_new, measurement_new, operator.le)
                cluster_new = Cluster({idx:sorted(list(component)) for idx, component in enumerate(nx.connected_components(G_new))})
                silhouette_score_new, silhouette_score_samples_new = cluster_new.silhouette(measurement_new)
                
                #print(f"{iter}\t{c}\t{node}\t{silhouette_score_new:.5f}\t{silhouette_score_tmp:.5f}\t{best_step['silhouette_score']:.5f}\t{silhouette_score_new > silhouette_score_tmp}\t{best_step_tmp['silhouette_score'] > best_step['silhouette_score']}")
                
                # best_step_tmp record the best step for each negative cluster
                if silhouette_score_new > silhouette_score_tmp:
                    best_step_tmp['graph'] = G_new
                    best_step_tmp['cluster'] = cluster_new
                    best_step_tmp['silhouette_score'] = silhouette_score_new
                    best_step_tmp['silhouette_score_samples'] = silhouette_score_samples_new
                    silhouette_score_tmp = silhouette_score_new
                
            # best_step choose the one to remove and update the graph
            if best_step_tmp['silhouette_score'] > best_step['silhouette_score']:
                best_step['graph'] = best_step_tmp['graph']
                best_step['cluster'] = best_step_tmp['cluster']
                best_step['silhouette_score'] = best_step_tmp['silhouette_score']
                best_step['silhouette_score_samples'] = best_step_tmp['silhouette_score_samples']

        result = Cluster({idx:sorted(list(component)) for idx, component in enumerate(nx.connected_components(best_step['graph']))})
        return result


    def optimize(self, sequences, measurement, method='ratio'):
        """
        Optimize the cluster

        Parameters
        ----------
        measurement : tuple
            List of measurement (seq1, seq2, measurement)
        method : str
            Optimization method
        
        Returns
        -------
        result : Cluster
            Optimized cluster
        """
        if method == 'ratio':
            return self._optimize_ratio(sequences, measurement)
        else:
            raise ValueError(f"Invalid optimization method: {method}")

           

