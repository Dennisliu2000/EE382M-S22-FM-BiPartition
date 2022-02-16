#######################################################################
# Implementation of FM partition
# You need to implement initialize() and partition_one_pass()
# All codes should be inside FM_Partition class
# Name: Dennis Liu
# UT EID: dl34437
#######################################################################

from typing import List, Tuple

import numpy as np

import math as m

from .p1_partition_base import FM_Partition_Base

__all__ = ["FM_Partition"]

class FM_Partition(FM_Partition_Base):
    def __init__(self) -> None:
        super().__init__()
    
    def update_gbuckets(self):                                  # Calculates gains for all nodes - technically only neighbors are needed, but this is what I have for now

        self.block0_g.clear()
        self.block1_g.clear()
        
        for n in range(len(self.partition[0])):                 # Computing gains for each node in block0 - moves each node in block0 to block1, calculates cut size diff, then moves them back
            if self.partition[0][n] in self.locked_nodes:       # If the current node we are looking at is in the locked nodes list, just ignore it
                continue
            else:
                self.partition[1].append(self.partition[0][n])      
                self.partition[0].pop(n)
                temp_cut_size = self.compute_cut_size(self.partition)
                temp_cut_size = self.cut_size - temp_cut_size
                self.partition[0].insert(n, self.partition[1][-1])
                self.partition[1].pop()
                if temp_cut_size in self.block0_g.keys():
                    self.block0_g.get(temp_cut_size).append(self.partition[0][n])
                    self.block0_g.get(temp_cut_size).sort()
                else:
                    self.block0_g[temp_cut_size] = [self.partition[0][n]]
        
        for n in range(len(self.partition[1])):                 # Computing gains for each node in block1
            if self.partition[1][n] in self.locked_nodes:
                continue
            else:
                self.partition[0].append(self.partition[1][n])
                self.partition[1].pop(n)
                temp_cut_size = self.compute_cut_size(self.partition)
                temp_cut_size = self.cut_size - temp_cut_size
                self.partition[1].insert(n, self.partition[0][-1])
                self.partition[0].pop()
                if temp_cut_size in self.block1_g.keys():
                    self.block1_g.get(temp_cut_size).append(self.partition[1][n])
                    self.block1_g.get(temp_cut_size).sort()
                else:
                    self.block1_g[temp_cut_size] = [self.partition[1][n]]
    
    ''' DEAD FUNCTION
    def update_gnodes(self):
        
        listk = []
        listk1 = []

        for key, lists in self.block0_g.items():
            for l_nodes in self.locked_nodes:
                if l_nodes in lists:
                    lists.remove(l_nodes)
                    if not lists:
                        listk.append(key)
                        
        for key in listk:
            self.block0_g.pop(key)
                    
        for key, lists in self.block1_g.items():
            for l_nodes in self.locked_nodes:
                if l_nodes in lists:
                    lists.remove(l_nodes)
                    if not lists:
                        listk1.append(key)
                    
        for key in listk1:
            self.block1_g.pop(key)
    '''
        
    def move_rtl(self):
        
        #x = sorted(self.block1_g.keys())                 # Get all gains in block1 and sort
        
        biggest_gains = self.block1_g.get(max(self.block1_g))        # Get the highest gain's list of nodes
        biggest_gains.sort()
        self.locked_nodes.append(biggest_gains[0])
        self.partition[0].append(biggest_gains[0])      # Add in the first node to the "left" (block0)
        self.partition[1].remove(biggest_gains[0])      # Remove the node from the "right" (block 1)
        
        
    def move_ltr(self):                                 # Same as move_rtl, but the blocks have been switched
        
        #x = sorted(self.block0_g.keys())
        
        biggest_gains = self.block0_g.get(max(self.block0_g))
        biggest_gains.sort()
        self.locked_nodes.append(biggest_gains[0])
        self.partition[1].append(biggest_gains[0])
        self.partition[0].remove(biggest_gains[0])
        
        
    def move(self):
        
        #x = sorted(self.block1_g.keys())  
        #y = self.block0_g.keys().sort()
        if self.n_nodes - len(self.locked_nodes) > 1:
            if max(self.block1_g) > max(self.block0_g):       # If the largest gain on block1 is greater than the largest gain in block0, perform a right-to-left move
                self.move_rtl()
            elif max(self.block0_g) > max(self.block1_g):     # Vice-versa
                self.move_ltr()
            else:                   # Both blocks have an equal maximum gain; find the first node via node number tie-breaker and then perform a rtl or ltr move
                a = self.block0_g.get(max(self.block0_g))
                b = self.block1_g.get(max(self.block1_g))
                a.sort()
                b.sort()
                if a[0] < b[0]:
                    self.move_rtl()
                else:
                    self.move_ltr()
        else:                                                   # There is only one node left to move
            if self.block0_g:
                self.move_ltr()
            else:
                self.move_rtl()
             
        
    def initialize(self):
        """Initialize necessary data structures before starting solving the problem
        """
        self.cut_size = 0
        self.cut_size_list = []
        self.best_sol = ()
        self.best_cut_size = 0
        self.locked_nodes = []
        self.block0 = []
        self.block1 = []
        self.r_ep = 0
        self.block0_g = {}
        self.block1_g = {}
        
        # TODO initial solutions: block 0 and block 1
        # To ensure a deterministic solution, use the following partition as the initial solution
        # sort the node names in alphabetical order
        # the first floor(N/2) nodes are in the first partition, The rest N-floor(N/2) nodes are in the second partition
        # a_0, a_1, ..., a_N/2-1 | a_N/2, a_N/2+1, ..., a_N-1, if N even
        # a_0, a_1, ..., a_(N-3)/2 | a_(N-1)/2, ..., a_N-1, if N odd
        # ...
        
        # Creates the initial partition
        if self.n_nodes % 2 == 0:
            for x in range(int(self.n_nodes / 2)):
                self.block0.append(x)
            for y in range(int(self.n_nodes / 2), self.n_nodes):
                self.block1.append(y)
            self.partition = (self.block0, self.block1)
        else:
            for x in range(m.floor(self.n_nodes / 2)):
                self.block0.append(x)
            for y in range(m.floor(self.n_nodes / 2), self.n_nodes):
                self.block1.append(y)
            self.partition = (self.block0, self.block1)                 # Partition is Tuple(List[int], List[int]) with int representation of nodes
            
        # Sets initial values after the first partition
        self.cut_size = self.compute_cut_size(self.partition)
        self.cut_size_list.append(self.cut_size)
        self.best_cut_size = self.cut_size
        self.r_ep = self.min_cut_ratio - self.min_cut_ratio_epsilon
        
        self.best_sol = ([self.node2node_name_map[n] for n in self.partition[0]], [self.node2node_name_map[n] for n in self.partition[1]])
        
        # Initialize gain buckets
        self.update_gbuckets()
        
        #print(self.partition)
        #print(self.block0_g)
        #print(self.block1_g)
        #print(self.cut_size)
            
        # TODO initialize any auxiliary data structure you need
        # e.g., node2net_map, cell gains, locked cells, etc.


    def partition_one_pass(self) -> Tuple[List[int], Tuple[List[str], List[str]], int]:
        """FM graph partition algorithm for one pass

        Return:
            cut_size_list (List[int]): contains the initial cut size and the cut size after each move
            best_sol (Tuple[List[str], List[str]]): The best partition solution is a tuple of two blocks.
                Each block is a list of node names. (Please use the original node names from the benchmark file.
                Hint: you might need to use node2node_name_map). If multiple solutions have the best cut size, return the first one.
            best_cut_size (int): The cut size corresponding to the best partition solution
        """
        # TODO implement your FM partition algorithm for one pass.
        # To make this method clean, you can extract subroutines as methods of this class
        # But do not override methods in the parent class
        # Please strictly follow the return type requirement.
        
        # (m.min(len(self.partition[0]), len(self.partition[1])) / self.n_nodes) < r_ep      
        while len(self.locked_nodes) < self.n_nodes:                        # As long as the locked-nodes list does not have all the nodes, we can iterate and swap
            print(self.partition)
            print(self.block0_g)
            print(self.block1_g)
            print(self.cut_size)
            self.rtl_min = min((len(self.partition[1]) - 1), (len(self.partition[0]) + 1))
            self.ltr_min = min((len(self.partition[0]) - 1), (len(self.partition[1]) + 1))
            if (self.rtl_min / self.n_nodes) < self.r_ep:   # If a right-to-left (block1 to block0) move would make the right partition too small, then execute a left-to-right move
                self.move_ltr()
            elif (self.ltr_min / self.n_nodes) < self.r_ep: # If a left-to-right move would make the left partition too small, then execute a right-to-left move
                self.move_rtl()
            else:                                                           # If both moves are acceptable, then find the best move out of both sides
                self.move()
                
            #self.update_gnodes()
            
            self.update_gbuckets()                                          # After a move, the gains are re-calculated
            self.cut_size = self.compute_cut_size(self.partition)           # Add the new cutsize to the history; check to see if we have a new best
            self.cut_size_list.append(self.cut_size)
            
            if self.cut_size < self.best_cut_size:
                self.best_cut_size = self.cut_size
                self.best_sol = ([self.node2node_name_map[n] for n in self.partition[0]], [self.node2node_name_map[n] for n in self.partition[1]])
                
            #print(self.partition)
            #print(self.block0_g)
            #print(self.block1_g)
            #print(self.cut_size)
                
        else: 
            return self.cut_size_list, self.best_sol, self.best_cut_size    # No moves available, first pass is over. Return gathered values.
        
