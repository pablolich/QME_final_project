import copy
from integration import *
import numpy as np

def index_mapping(old_ind, del_ind):
    '''
    Given lists of indices of certain positions and deletions on a vector, 
    determine the new indices of positions once deletions are removed.
    Note that the intersection between old_ind and del_ind must be the empty
    set, and also that their union need not span the full length of the vector.

    Example:

        vector = np.array([1, 2, 3, 4, 5])
        old_index = [0, 3]
        del_index = [1, 4]
        new_index = index_mapping(old_index, del_index)
        print(vector[old_index])
        new_vector = np.delete(vector, del_index)
        print(new_vector[new_index])
        #the two print statements yield the same output
    '''
    return [i - sum([j < i for j in del_ind]) for i in old_ind]


class Community:
    def __init__(self, n, A, r, model):
        self.n = n #abundance vector
        self.A = A #competition matrix
        self.r = r #growth rate
        self.model = model #model on which the instance of the class is based
        self.presence = np.zeros(len(n), dtype = bool)
        ind_extant = np.where(n > 0)[0]
        self.presence[ind_extant] = True #indices of present species
        self.richness = len(ind_extant) #richness

    def remove_spp(self, remove_ind):
        '''
        remove all species in vector 'remove_ind' from community
        '''
        if self.model.__name__ == 'GLV':
            #create a deep copy of comm to keep original unmodified
            new_comm = copy.deepcopy(self)
            #remove row and column indices 'remove_ind' from A
            del_row = np.delete(new_comm.A, remove_ind, axis=0)
            del_row_col = np.delete(del_row, remove_ind, axis=1)
            new_comm.A  = del_row_col
            #remove element from abundance and growth rate vectors
            new_comm.n = np.delete(new_comm.n, remove_ind)
            new_comm.r = np.delete(new_comm.r, remove_ind)
            #update presence vector
            new_comm.presence[remove_ind] = False
            #get number of species actually removed (i.e., only those whose 
            #abundance was different than 0)
            n_rem = sum(self.n[remove_ind]>0)
            #reduce richness accordingly
            new_comm.richness -= n_rem
        else:
            raise ValueError('unknown model name')
        return new_comm

    def add_spp(self, add_ind, **kwargs):
        '''
        add all the species in 'add_ind' which details are in **kwargs
        '''
        if self.model.__name__ == 'GLV':
            #create a deep copy of comm to keep original unmodified
            new_comm = copy.deepcopy(self)
            #switch to ones in the indices of introduced species
            new_comm.presence[add_ind] = True
            mask = new_comm.presence == True
            add_row = kwargs['row'][mask]
            add_col = kwargs['col'][mask]
            #map old index vector into new index vector
            new_add = index_mapping(add_ind, 
                                    np.where(new_comm.presence==False)[0])
            #update richness
            new_comm.richness += len(new_add)
            #delete diagonal element to adhere to previous dimensions of A
            add_row_d = np.delete(add_row, new_add)
            #add rows and columns at the end of matrix A
            new_comm.A = np.insert(new_comm.A, new_add, add_row_d, axis = 0)
            new_comm.A = np.insert(new_comm.A, new_add, 
                                   add_col.reshape(new_comm.richness, 
                                                   len(new_add)), axis = 1)
            #add element to growth rate
            new_comm.r = np.insert(new_comm.r, new_add, kwargs['r'])
            #update abundances
            new_comm.n = np.insert(new_comm.n, new_add, kwargs['x'])
        else:
            raise ValueError('unknown model name')
        return new_comm

    def is_subcomm(self, presence):
        '''
        determine if the presence/absence binary vector is a subset of the 
        community
        '''
        #CHECK IF THE VECTOR PRESENCE HAS TO BE BINARY OR IT CAN BE BOOLEAN
        set1 = set(np.where(self.presence == True)[0])
        set2 = set(np.where(presence == 1)[0])
        if set1 == set2:
            return False
        else: 
            return set1.issubset(set2)

    def assembly(self, tol=1e-9, delete_history=False):
        if self.model.__name__ == 'GLV':
            #integrate using lemke-howson algorithm
            n_eq = lemke_howson_wrapper(-self.A, self.r)
            if np.all(n_eq == 0):
                #lemke-howson got stuck, integrate using differential equations
                print('integrating dynamics the hard way')
                sol = prune_community(self.model, self.n, tol,
                                      args=(self.A, self.r), 
                                      events=single_extinction)
                n_eq = sol.y[:,-1]
            #set to 0 extinct species
            ind_ext = np.where(n_eq < tol)[0]
            if any(ind_ext):
                n_eq[ind_ext] = 0
            self.n = n_eq
            self.presence[ind_ext] = False
            self.richness -= len(ind_ext) #update richness 
        else:
            print("haven't coded up other type of models yet")
            raise ValueError
        return self

    def delete_history(self):
        '''
        Delete history of assemlby, that is remove zeroed species, as well as
        absences from the presence vector
        '''
        #remove extinct species
        rem_ind = np.where(self.presence == 0)[0]
        comm = self.remove_spp(rem_ind)
        #remove from presence vector
        comm.presence = self.presence[self.presence]
        return comm

class Environment:
    def __init__(self, r):
        self.r = r #supply rate of each resource
