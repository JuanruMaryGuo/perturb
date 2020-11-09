import torch
from torch.utils.data import Dataset
from numpy.random import choice as npc
import numpy as np
import random
import helpfuntions


class perturbdataloader(Dataset):

    def __init__(self, dataset, batchsz, n_ways = 10, k_shots = 5, k_query = 15, shuffle = True , plus = 0):
        super(perturbdataloader, self).__init__()
        self.dataset = dataset
        self.batchsz = batchsz
        self.n_ways = n_ways
        self.k_shots= k_shots
        self.k_query = k_query
        self.shuffle = shuffle
        self.plus = plus

    def __len__(self):
        return self.batchsz

    def __getitem__(self, index):
        select = random.sample(range(self.n_ways), self.n_ways)

        self.datas = helpfuntions.make_datas(self.dataset, select,self.plus)

        inputs_support, inputs_query, target_support, target_query = helpfuntions.sample_once(self.datas,
                                                                                              self.k_shots,
                                                                                              self.k_query,
                                                                                              shuffle_ornot =self.shuffle,
                                                                                              plus = 0)

        return torch.from_numpy(inputs_support).float(), torch.from_numpy(target_support).long(), \
               torch.from_numpy(inputs_query).float(), torch.from_numpy(target_query).long()

