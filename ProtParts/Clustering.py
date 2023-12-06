import networkx as nx
import operator
from sklearn.metrics import silhouette_score, silhouette_samples
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
 
    
    def evaluate_matrix(self, method=None, matrix=None):
        """
        Evaluate the cluster based on similarity/distance matrix

        Parameters
        ----------
        method : str
            Evaluation method
        matrix : np.array
            Similarity/Distance matrix
        
        Returns
        -------
        metric : float
            Evaluation result
        """

        if method == 'sihouette':
            return self._silhouette(matrix)

    
    def _silhouette(self, matrix):
        """
        Silhouette score

        Parameters
        ----------
        matrix : np.array
            Similarity/Distance matrix
        
        Returns
        -------
        metric : float
            Silhouette score
        """
        pass




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
            clusters = self._graph(sequences, measurement, op)
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
        measurement = list(filter(lambda x:(x[0] != x[1]) and (op(x[2], self.threshold)), measurement))
        
        # add edges
        G.add_weighted_edges_from(measurement)
        
        result = {idx:sorted(list(component)) for idx, component in enumerate(nx.connected_components(G))}

        return Cluster(result)
    

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

                

