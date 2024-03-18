import random
# import time

# class Partition:
#     """
#     Partition class
#     """
#     def __init__(self, partitions=None, num_partitions=None):
#         """
#         Parameters
#         ----------
#         partition : dict
#             Partitioned clusters, {0: {1: [seq1, seq2, ...], 2: [seq1, seq2, ...], ...}, ...}
#         """
#         if partitions:
#             self.partitions = partitions
#         elif num_partitions:
#             self.partitions = {i:dict() for i in range(num_partitions)}
#         else:
#             raise ValueError('Either partitions or num_partitions must be specified')


#     def __len__(self):
#         """
#         Returns
#         -------
#         length : int
#             Number of partitions
#         """
#         return len(self.partitions)
    

#     def __getitem__(self, key):
#         """
#         Parameters
#         ----------
#         key : int
#             Partition ID
        
#         Returns
#         -------
#         partition : dict
#             Partitioned clusters
#         """
#         return self.partitions[key]
    

#     def __setitem__(self, key, value):
#         """
#         Parameters
#         ----------
#         key : int
#             Partition ID
#         value : dict
#             Partitioned clusters
#         """
#         self.partitions[key] = value

    
#     def items(self):
#         """
#         Returns
#         -------
#         items : list
#             List of partition items
#         """
#         return self.partitions.items()


#     def num_cluster(self, by='sum'):
#         """
#         Parameters
#         ----------
#         by : str/int
#             'sum': sum of clusters in all partitions
#             'max': max of clusters among partitions
#             'min': min of clusters among partitions
#             str/int: number of clusters in the partition with the given partition ID
        
#         Returns
#         -------
#         size : int
#             Number of clusters in the partition
#         """
#         size_list = [len(p) for p in self.partitions.values()]
#         if by == 'sum':
#             size = sum(size_list)
#         elif by == 'max':
#             size = max(size_list)
#         elif by == 'min':
#             size = min(size_list)
#         elif isinstance(by, int) or isinstance(by, str):
#             size = len(self.partitions[by])
        
#         return size
    

#     def num_data(self, by='sum'):
#         """
#         Parameters
#         ----------
#         by : str/int
#             'sum': sum of data in all partitions
#             'max': max of data among partitions
#             'min': min of data among partitions
#             'all': list of data in all partitions
#             str/int: number of data in the partition with the given partition ID
        
#         Returns
#         -------
#         size : int/list
#             Number of data in the partition/List of data in all partitions
#         """
#         size_list = [sum(len(c) for c in p.values()) for p in self.partitions.values()]
#         if by == 'sum':
#             size = sum(size_list)
#         elif by == 'max':
#             size = max(size_list)
#         elif by == 'min':
#             size = min(size_list)
#         elif isinstance(by, int) or isinstance(by, str):
#             size = sum(len(c) for c in self.partitions[by].values())
#         elif by == 'all':
#             size = size_list
        
#         return size






class Partitioning:
    """
    Partitioning class
    """
    def __init__(self, num_partitions, num_sequences, method='random'):
        """
        Parameters
        ----------
        num_partitions : int
            Number of partitions
        method : str
            Partitioning method
        """
        self.num_partitions = num_partitions
        self.num_sequences = num_sequences
        self.method = method

    
    def random_partitioning(self, clusters, random_seed=0):
        """
        Random partitioning of clusters

        Parameters
        ----------
        cluster : Cluster
            Cluster object
        random_seed : int
            Random seed
        
        Returns
        -------
        partitions : list
            Dict of partitions
        """
        random.seed(random_seed)

        # print([len(v) for k, v in clusters.items()])
        
        #partitions = Partition(num_partitions=self.num_partitions)
        partitions = {i:dict() for i in range(self.num_partitions)}
        num_sample = clusters.num_data(by='sum')
        if num_sample != self.num_sequences:
            raise ValueError(f"Number of sequences in clusters {num_sample} is different from the number of sequences {self.num_sequences}")
        partitions_alloc_size = self._even_split()
        done_list = []
        # print(partitions_alloc_size)  

        # sort clusters by size
        clusters_sorted = clusters.sort(reverse=True)

        if clusters.num_data(by='max') > max(partitions_alloc_size.values()):
            raise ValueError(f"Cluster size {clusters.num_data(by='max')} is larger than maximum partition size {max(partitions_alloc_size)}")
        

        for c_id, c in clusters_sorted.items():
            merged = False
            # print(c_id, len(c))
            # print partition size
            left_alloc_size = {key: sum(len(sublist) for sublist in value.values()) for key, value in partitions.items()}
            # print(left_alloc_size)
            while not merged:
                #partition_id = random.randint(0, self.num_partition-1)
                # candicate idx is the key of partition which has a number less than the corresponding partition size
                candidate_idx = [i for i in range(self.num_partitions) if i not in done_list]
                partition_id = random.choice(candidate_idx)
                #partition_id_key = f"Partition_{partition_id}"
                # print('\t', partition_id, left_alloc_size[partition_id], len(c))
                if left_alloc_size[partition_id] + len(c) <= partitions_alloc_size[partition_id]:
                    # partition[partition_id].append(clusters_id)
                    partitions[partition_id][c_id] = c
                    merged = True
                    # print(merged, '\n')
                    if left_alloc_size[partition_id] + len(c) == partitions_alloc_size[partition_id]:
                        done_list.append(partition_id)
                # time.sleep(0.1)
                

        return partitions

    
    def _even_split(self):
        """
        Evenly split clusters into partitions

        Returns
        -------
        par_size : dict
            Dict of partition size
        """
        par_size = [int(self.num_sequences/self.num_partitions)]*self.num_partitions
        for i in range(self.num_sequences%self.num_partitions):
            par_size[i] += 1
        return {i:size for i, size in enumerate(par_size)}
    

    def partition_size(self):
        """
        Calculate partition size

        Returns
        -------
        partition_size : dict
            Dict of number of sequences in each partition
        """
        return self._even_split()






