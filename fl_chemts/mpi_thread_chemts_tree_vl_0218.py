from subprocess import Popen, PIPE
from math import *
import random
import numpy as np
import random as pr
from copy import deepcopy
#from types import IntType, ListType, TupleType, StringTypes
import itertools
import time
import math
import os
import shutil

import tensorflow as tf

import argparse
import subprocess
from load_model import loaded_model
from keras.preprocessing import sequence
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors
from rdkit.Chem import MolFromSmiles, MolToSmiles
from rdkit.Chem import Crippen
import sys
#from make_smile import zinc_data_with_bracket_original, zinc_processed_with_bracket
from threading import Thread, Lock, RLock
import threading
from Queue import Queue
from mpi4py import MPI
from RDKitText import tansfersdf
from SDF2GauInput import GauTDDFT_ForDFT
from GaussianRunPack import GaussianDFTRun

import sascorer
import pickle
import gzip
import networkx as nx
from rdkit.Chem import rdmolops


smiles_max_len = 81 # zinc dataset
state_length = 64


class chemical:

    def __init__(self):

        self.position=['&']
    def Clone(self):

        st = chemical()
        st.position= self.position[:]
        return st

    def SelectPosition(self,m):

        self.position.append(m)

    def Getatom(self):
        return [i for i in range(self.num_atom)]

class Node:

    def __init__(self, position = None, parent = None, state = None, nodelock=threading.Lock()):
        self.position = position
        self.parentNode = parent
        self.childNodes = []
        self.child=None
        self.wins = 0
        self.re_max = 0
        self.visits = 0
        self.depth=0
        self.expanded=[]
        self.nodeadded=[]
        self.random_node=[]
        self.all_posible=[]
        self.generate_smile=[]
        self.node_index=[]
        self.valid_smile=[]
        self.new_compound=[]
        self.nodelock=nodelock
        self.ucb=[]
        self.core_id=[]
        self.virtual_loss=0
        self.num_thread_visited=0
        self.all_probs=[]


    def Selectnode(self, ts_strategy, search_parameter, alpha):
        #self.nodelock.acquire()
 
        ucb=[]
        ntv_list=[]
        base_list=[]
        bias_list=[]
        max_list=[]
        #print "current node's virtual_loss:",self.num_thread_visited,self.virtual_loss
        for i in range(len(self.childNodes)):
            #print "current node's childrens' virtual_loss:",self.childNodes[i].num_thread_visited,self.childNodes[i].virtual_loss
            C = search_parameter
            cNodei = self.childNodes[i]
            if ts_strategy == 'uct': 
                ucb.append(alpha*(cNodei.wins)/(0.0001+cNodei.visits+cNodei.num_thread_visited)+
                           (1-alpha)*cNodei.re_max/(1+cNodei.num_thread_visited)+
                           C*sqrt(2*log(self.visits+self.num_thread_visited)/(0.0001+cNodei.visits+cNodei.num_thread_visited)))
            elif ts_strategy == 'puct':
                prob=self.all_probs[i]
                ucb.append(alpha*(cNodei.wins)/(0.001+cNodei.visits+cNodei.num_thread_visited)+
                           (1-alpha)*cNodei.re_max/(1+cNodei.num_thread_visited)+
                           C*(np.tanh(2*prob-1)+1)/2*sqrt((self.visits+self.num_thread_visited))/(1+cNodei.visits+cNodei.num_thread_visited))
            ntv_list.append(cNodei.num_thread_visited)
            base_list.append(alpha*(cNodei.wins)/(0.001+cNodei.visits+cNodei.num_thread_visited)+(1-alpha)*cNodei.re_max/(1+cNodei.num_thread_visited))
            bias_list.append(ucb[-1] - base_list[-1])
            max_list.append(cNodei.re_max)
        #print 'ucb score list', ucb
        #print 'ntv_list', ntv_list, 'cNodei.num_thread_visited', cNodei.num_thread_visited, 'total', np.sum(ntv_list)
        #print 'base_list', base_list
        #print 'bias_list', bias_list
        #print 'max_list', max_list

        m = np.amax(ucb)
        indices = np.nonzero(ucb == m)[0]
        ind=pr.choice(indices)
        s=self.childNodes[ind]
        #print "which thread's ucb:",threading.currentThread().getName()
        #print ucb
        #self.nodelock.release()
        return s

    def Addnode(self, m):

        #n = Node(position = m, parent = self, state = s)
        self.nodeadded.remove(m)
        n = Node(position = m, parent = self)
        self.childNodes.append(n)
        return n



    def Update(self, result, add_vis_count = 1):
        #self.nodelock.acquire()
        #print "update visits:",self.visits
        self.visits += add_vis_count
        self.wins += result
        if self.re_max < result:
            self.re_max = result
        #self.nodelock.release()

    def delete_virtual_loss(self):
        #self.num_thread_visited=0
        self.num_thread_visited += -1
        self.virtual_loss=0

    def expanded_node1(self, model, state, val):


        all_nodes=[]

        end="\n"
        position=[]
        position.extend(state)
        total_generated=[]
        new_compound=[]
        get_int_old=[]
        for j in range(len(position)):
            get_int_old.append(val.index(position[j]))

        get_int=get_int_old
        x=np.reshape(get_int,(1,len(get_int)))
        #x_pad= sequence.pad_sequences(x, maxlen=42, dtype='int32', padding='post', truncating='pre', value=0.) #original
        x_pad= sequence.pad_sequences(x, maxlen=82, dtype='int32', padding='post', truncating='pre', value=0.)  #zinc 250,000 
        ex_time=time.time()

        for i in range(1):
            global graph
            with graph.as_default():
                predictions=model.predict(x_pad)
                #print "shape of RNN",predictions.shape
                preds=np.asarray(predictions[0][len(get_int)-1]).astype('float64')
                preds = np.log(preds) / 1.0
                preds = np.exp(preds) / np.sum(np.exp(preds))
                #next_probas = np.random.multinomial(1, preds, 1)
                #print('preds', preds)
                next_probas=np.argsort(preds)[-5:]
                next_probas=list(next_probas)
                #print('next_probas', next_probas)
		#next_int=np.argmax(next_probas)
                #get_int.append(next_int)
                #all_nodes.append(next_int)

        #all_nodes=list(set(all_nodes))
        if 0 in next_probas:
            next_probas.remove(0)
        all_nodes=next_probas
	#print('all_nodes', all_nodes)

        self.expanded=all_nodes
        #print self.expanded
        exfi_time=time.time()-ex_time
	#print exfi_time


    def expanded_node(self, model,state,val):
        all_nodes=[]

        end="\n"
        position=[]
        position.extend(state)
        total_generated=[]
        new_compound=[]
        get_int_old=[]
        for j in range(len(position)):
            get_int_old.append(val.index(position[j]))

        get_int=get_int_old
        x=np.reshape(get_int,(1,len(get_int)))
        x_pad= sequence.pad_sequences(x, maxlen=smiles_max_len, dtype='int32',
            padding='post', truncating='pre', value=0.)
	    #ex_time=time.time()
        for i in range(60):
            global graph
            with graph.as_default():
                predictions=model.predict(x_pad)
                #print "shape of RNN",predictions.shape
                preds=np.asarray(predictions[0][len(get_int)-1]).astype('float64')
                preds = np.log(preds) / 1.0
                preds = np.exp(preds) / np.sum(np.exp(preds))
                #print('preds', preds)
                next_probas = np.random.multinomial(1, preds, 1)
                next_int=np.argmax(next_probas)
                #print('next_int', next_int)
                #get_int.append(next_int)
                all_nodes.append(next_int)

        all_nodes=list(set(all_nodes))
    
        #print('all_nodes', all_nodes)
        self.expanded=all_nodes
        #print self.expanded
	#exfi_time=time.time()-ex_time
	#print exfi_time


    def expanded_node_puct(self, model,state,val):
        all_nodes=[]

        end="\n"
        position=[]
        position.extend(state)
        total_generated=[]
        new_compound=[]
        get_int_old=[]
        for j in range(len(position)):
            get_int_old.append(val.index(position[j]))

        get_int=get_int_old
        x=np.reshape(get_int,(1,len(get_int)))
        #x_pad= sequence.pad_sequences(x, maxlen=42, dtype='int32',
        #    padding='post', truncating='pre', value=0.)
        x_pad= sequence.pad_sequences(x, maxlen=smiles_max_len, dtype='int32',padding='post', truncating='pre', value=0.)
        #ex_time=time.time()                                                                                                               
        for i in range(1):
            global graph
            with graph.as_default():
                predictions=model.predict(x_pad)
                #print "shape of RNN",predictions.shape                                                                                        
                preds=np.asarray(predictions[0][len(get_int)-1]).astype('float64')
                preds = np.log(preds) / 1.0
                preds = np.exp(preds) / np.sum(np.exp(preds))
                #print('preds', preds, len(preds))
                next_probas = np.random.multinomial(1, preds, 1)
                next_int=np.argmax(next_probas)
                #print('next_int', next_int)
                #get_int.append(next_int)                                                                                                     
                #all_nodes.append(next_int)
        
        ordered_preds = np.sort(preds)[::-1]
        ordered_index = np.argsort(preds)[::-1]
        #print('ordered_preds', ordered_preds, 'ordered_index', ordered_index)
        cut_index = 0
        p_sum = 0
        for i in range(len(ordered_preds)):
            p_sum += ordered_preds[i]
            #print(i, p_sum)
            if p_sum > 0.99:
                cut_index = i+1
                break
        #all_nodes=list(set(all_nodes))
        all_nodes = ordered_index[:cut_index]
        all_probs = ordered_preds[:cut_index]
        #print('all_nodes', all_nodes, 'all_probs', all_probs)
        self.expanded=all_nodes
        self.all_probs=all_probs



    def node_to_add(self, all_nodes,val):
        added_nodes=[]
        #print('val',val)
        for i in range(len(all_nodes)):
            #print('val[all_nodes[i]]',val[all_nodes[i]],)
            added_nodes.append(val[all_nodes[i]])

        self.nodeadded=added_nodes

        #print "childNodes of current node:", self.nodeadded

    def random_node_to_add(self, all_nodes,val):
        added_nodes=[]
        for i in range(len(all_nodes)):
            added_nodes.append(val[all_nodes[i]])

        self.random_node=added_nodes





        #print "node.nodeadded:",self.nodeadded





"""Define some functions used for RNN"""



def chem_kn_simulation(model,state,val,added_nodes):
    all_posible=[]

    end="\n"

    position=[]
    position.extend(state)
    position.append(added_nodes)
    total_generated=[]
    new_compound=[]
    get_int_old=[]
    for j in range(len(position)):
        get_int_old.append(val.index(position[j]))

    get_int=get_int_old

    x=np.reshape(get_int,(1,len(get_int)))
    x_pad= sequence.pad_sequences(x, maxlen=smiles_max_len, dtype='int32',
        padding='post', truncating='pre', value=0.)
    while not get_int[-1] == val.index(end):
        predictions=model.predict(x_pad)
        #print "shape of RNN",predictions.shape
        preds=np.asarray(predictions[0][len(get_int)-1]).astype('float64')
        preds = np.log(preds) / 1.0
        preds = np.exp(preds) / np.sum(np.exp(preds))
        next_probas = np.random.multinomial(1, preds, 1)
        next_int=np.argmax(next_probas)
        a=predictions[0][len(get_int)-1]
        next_int_test=sorted(range(len(a)), key=lambda i: a[i])[-10:]
        get_int.append(next_int)
        x=np.reshape(get_int,(1,len(get_int)))
        x_pad = sequence.pad_sequences(x, maxlen=smiles_max_len, dtype='int32',
            padding='post', truncating='pre', value=0.)
        if len(get_int)>state_length:
            break
    total_generated.append(get_int)
    all_posible.extend(total_generated)


    return all_posible




def predict_smile(all_posible,val):
    new_compound=[]
    for i in range(len(all_posible)):
        total_generated=all_posible[i]

        generate_smile=[]

        for j in range(len(total_generated)-1):
            generate_smile.append(val[total_generated[j]])
        generate_smile.remove("&")
        new_compound.append(generate_smile)

    return new_compound


def make_input_smile(generate_smile):
    new_compound=[]
    for i in range(len(generate_smile)):
        middle=[]
        for j in range(len(generate_smile[i])):
            middle.append(generate_smile[i][j])
        com=''.join(middle)
        new_compound.append(com)
    return new_compound




def ChemTS_run(rootnode,result_queue,lock,chem_model,ts_strategy,search_parameter,num_simulations, gau_parallel,simulation_time, output_file,alpha,objective,num_rollout,charge_check,SA_score_check):
    """----------------------------------------------------------------------"""
    """----------------------------------------------------------------------"""
    global maxnum
    global gau_file_index
    global ind_mol
    start_time=time.time()
    #while time.time()-start_time<3600*24:
    while time.time()-start_time<simulation_time:
        node = rootnode
        state=['&']
        """selection step"""
        node_pool=[]
        lock.acquire()

        #print 'node.expanded', node.expanded, 'node.nodeadded', node.nodeadded, 'len(node.childNodes)', len(node.childNodes), len(node.expanded)
        while len(node.expanded)>0 and node.nodeadded==[] and len(node.childNodes)==len(node.expanded):
            #node.num_thread_visited+=1
            #node.virtual_loss+=0
            node = node.Selectnode(ts_strategy, search_parameter, alpha)
            state.append(node.position)
            #print 'node.expanded', node.expanded, 'node.nodeadded',node.nodeadded,'len(node.childNodes)',len(node.childNodes), len(node.expanded)
        depth.append(len(state))
	#lock.release()

        """this if condition makes sure the tree not exceed the maximum depth"""
        if len(state)>state_length:
            #lock.acquire()
            #dest_core=random.choice(free_core_id)
            #free_core_id.remove(dest_core)
            #comm.send([state,m], dest=dest_core, tag=START)
            #lock.release()
            #data = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
            #lock.acquire()
            #free_core_id.append(data[2])
            #lock.release()

            #tag = status.Get_tag()
            #if tag == DONE:
                #lock.acquire()
                #all_compounds.append(data[1])
                #lock.release()
            re=-10


            #lock.acquire()
            while node != None:
                node.Update(re)
                #node.delete_virtual_loss()
                node = node.parentNode
            lock.release()
        else:
            """expansion step"""
            #lock.acquire()
            #print('node.expanded', node.expanded, 'node.nodeadded', node.nodeadded)
            m = None
            if node.expanded==[]:
                #maxnum+=1
                #ind_mol+=1
                if ts_strategy == 'uct':
                    node.expanded_node(chem_model,state,val)
                elif ts_strategy == 'puct':
                    node.expanded_node_puct(chem_model,state,val)
                node.node_to_add(node.expanded,val)
                node.random_node_to_add(node.expanded,val)
                
                if node.nodeadded!=[]:
                    
                    m=node.nodeadded[0]
                    #node=node.Addnode(m)
            else:
                if node.nodeadded!=[]:
                    m=node.nodeadded[0]
                    #node=node.Addnode(m)

            if m == None:
                m = val[random.choice(node.expanded)]
                print('randomly selected')
            else:
                if m != '\n':
                    node = node.Addnode(m)
                else:
                    #print('node.nodeadded', node.nodeadded)
                    node.nodeadded.remove(m)
                    #print('node.nodeadded', node.nodeadded)
                    #print('\\n was tried to added... skip')
                    lock.release()   
                    continue
            #print "m is:",m

            lock.release()
	   
	    """simulation step"""
            for ro in range(num_rollout):
                                
                lock.acquire()
                """add virtual loss"""
                node_tmp = node
                while node_tmp != None:
                    #print "node.parentNode:",node.parentNode                                                                                                
                    #node.Update(re)
                    #node.delete_virtual_loss()
                    node_tmp.num_thread_visited+=1
                    node_tmp = node_tmp.parentNode
                print 'rootnode.num_thread_visited', rootnode.num_thread_visited
                lock.release()

                lock.acquire()
                maxnum+=1
                ind_mol+=1
                #print('ind_mol', ind_mol)
                #lock.release()
                #"""simulation step"""
                #lock.acquire()
                print('free_core_id_prev', len(free_core_id),'use_core_id', len(use_core_id))

                dest_core=random.choice(free_core_id)
                use_core_id.append(dest_core)
                free_core_id.remove(dest_core)
                print('dest_core', dest_core)
                #lock.release()
                try:
                    comm.send([state,m,ind_mol], dest=dest_core, tag=START)
                    lock.release()
                except:
                    print('comm.send failed')
                    free_core_id.append(dest_core)
                    use_core_id.remove(dest_core)
                    lock.acquire()
                    """backpropation step"""
                    while node!= None:
                        #print "node.parentNode:",node.parentNode                                                                                                              
                        node.Update(0, add_vis_count = 0)
                        node.delete_virtual_loss()
                        node = node.parentNode
                    lock.release()
                    
                    continue
                #lock.release()

                try:
                    #data = [uv, compound, rand, uv intensity, deen, gap]
                    #data = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
                    data = comm.recv(source=dest_core, tag=MPI.ANY_TAG, status=status)   

                    lock.acquire()
                    free_core_id.append(data[2])
                    use_core_id.remove(data[2])
                    print('data[2]', data[2], 'dest_core', dest_core)
                    lock.release()
                except:
                    print('comm.recv failed.')
                    lock.acquire()
                    #data = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
                    #print('except data', data)
                    free_core_id.append(dest_core)
                    use_core_id.remove(dest_core)
                    
                    #data = [-1000, '', 0, 0, 0, 0, 0]
                    """backpropation step"""
                    while node!= None:
                        #print "node.parentNode:",node.parentNode                       
                        node.Update(0, add_vis_count = 0)
                        node.delete_virtual_loss()
                        node = node.parentNode
                    lock.release()
                    
                    continue
                #lock.release()
                print('free_core_id', free_core_id)
			#with open('/home/yang/DP-ChemTS/csvresult.csv','wb') as file:
			 #   for line in text:
			#	    file.write(line)
			#		file.write('\n')
                #lock.acquire()
			#result_test.append[data]
                #free_core_id.append(data[2])
                #with open('/home/yang/csvfile.csv','wb') as file:
                #    for line in result_test:
                #        file.write(str(line))
                #        file.write('\n')
                #lock.release()
                re = 0
                tag = status.Get_tag()
                if tag == DONE:
                    lock.acquire()
                    all_compounds.append(data[1])
                    lock.release()
                    if data[0]!=-1000 and data[3] >= 0:
                        #if data[3] > 0.01:
                        #    re=(0.01*data[0])/(1.0+abs(0.01*data[0]))
                        #else:
                        if objective == 'WL_IT':
                            if data[0] < 700:
                                re = (np.tanh(0.003*(data[0]-400)) + 1)/2
                            else:
                                wl_re = (np.tanh(0.003*(data[0]-400)) + 1)/2
                                intensity_re = (np.tanh((np.log10(data[3]+0.00000001)-np.log10(0.01)))+1)/2
                                w_wl = 0.8
                                w_intensity = 0.2
                                re = w_wl *wl_re+ w_intensity *intensity_re
                        elif objective == 'HL':
                            #HOMO/LUMO                                                                                                                            
                            re = 1 - data[5]/10.
                        elif objective == 'WL':
                            re = (np.tanh(0.003*(data[0]-400)) + 1)/2
                        elif objective == 'WLFL':
                            fl_target = 1200
                            w_fl = 0.8
                            w_fl_intensity = 0.2
                            wl_target = 700
                            w_wl = 0.8
                            w_wl_intensity = 0.2
                            
                            wl = data[0]
                            wl_intensity = data[3]
                            fl = data[13]
                            fl_intensity = data[14]
                            

                            sigma = 150
                            wl_re = np.exp(-(wl-wl_target)**2/(2*sigma**2))
                            wl_intensity_re = (np.tanh((np.log10(wl_intensity+0.00000001)-np.log10(0.01)))+1)/2
                            re_wl = w_wl *wl_re+ w_wl_intensity *wl_intensity_re
                            
                            if fl == -1000:
                                re_fl = 0
                            else:
                                sigma = 100
                                fl_re = np.exp(-(fl - fl_target)**2/(2*sigma**2))
                                fl_intensity_re = (np.tanh((np.log10(np.abs(fl_intensity)+0.00000001)-np.log10(0.01)))+1)/2
                                re_fl = w_fl * fl_re+ w_fl_intensity * fl_intensity_re
                            
                            re = (re_wl+re_fl)/2


                        if SA_score_check:
                            data[11]
                            re = re*((-np.tanh(data[11]-4)+1)/2)
                        #re= (0.01*data[0])/(1.0+abs(0.01*data[0]))    

                        #For penality of duplication 
                        if data[1] in wave_compounds:
                            re = -1

                        lock.acquire()
                        wave_compounds.append(data[1])
                        wave.append(data[0])
                        deen_list.append(data[4])
                        uv_intensity_list.append(data[3])
                        gap_list.append(data[5])
                        wl_list_list.append(data[6])
                        intensity_list_list.append(data[7])
                        reward_list.append(re)
                        index_list.append(data[8])
                        mol_weight_list.append(data[9])
                        logP_list.append(data[10])
                        SA_score_list.append(data[11])
                        depth_list.append(data[12])
                        s1_wavelength_list.append(data[13])
                        s1_strength_list.append(data[14])
                        s1_wavelength_list_list.append(data[15])
                        s1_strength_list_list.append(data[16])

                        with open('csvcom_.csv','wb') as file: 
                        #with open('/home/yang/csvcom.csv','wb') as file:
                            for line1 in wave_compounds:
                                file.write(str(line1))
                                file.write('\n')
                        with open('csvwave_.csv','wb') as file:
                        #with open('/home/yang/csvwave.csv','wb') as file:
                            for line2 in wave:
                                file.write(str(line2))
                                file.write('\n')
                    
                        with open(output_file,'wb') as file:
                            file.write('#Search strategy, '+ts_strategy+', search_parameter, '+str(search_parameter)+', alpha, '+str(alpha)+', objective,'+str(objective)+', parallel simulations, '+str(num_simulations)+', gaussian parallel, '+str(gau_parallel)+', simulation_time (h), '+str(simulation_time/3600)+', num_rollout, '+str(num_rollout)+'charge_check, '+str(charge_check)+'SA_score_check, '+str(SA_score_check)+'\n')
                            file.write('#Compound, index, wavelength, uv_intensity, s1_wavelength, s1_strength, reward,  deen, gap, mol_weight, logP, SA_score, TS depth,  wavelength_list, intensity_list, s1_wl_list, s1_strength_list \n')
                            for i in range(len(wave_compounds)):
                                file.write(str(wave_compounds[i])+', ')
                                file.write(str(index_list[i])+', ')
                                file.write(str(wave[i])+', ')
                                file.write(str(uv_intensity_list[i])+', ')
                                file.write(str(s1_wavelength_list[i])+', ')
                                file.write(str(s1_strength_list[i])+', ')
                                file.write(str(reward_list[i])+', ')
                                file.write(str(deen_list[i])+', ')
                                file.write(str(gap_list[i])+', ')
                                file.write(str(mol_weight_list[i])+', ')
                                file.write(str(logP_list[i])+', ')
                                file.write(str(SA_score_list[i])+', ')
                                file.write(str(depth_list[i])+', ')
                                for wl_i in wl_list_list[i]:
                                    file.write(str(wl_i)+', ')
                                for int_i in intensity_list_list[i]:
                                    file.write(str(int_i)+', ')
                                for wl_i in s1_wavelength_list_list[i]:
                                    file.write(str(wl_i)+', ')
                                for int_i in s1_strength_list_list[i]:
                                    file.write(str(int_i)+', ')
                                file.write('\n')

                            
                            
                        lock.release()
                    if data[0]==-1000:
                        #re=-1
                        re=0
                    if data[3]<0:
                        #re=-1
                        re=0
                    #if m=='\n':
                    #    re=-10000

                lock.acquire()
                if re == None:
                    re = 0
                """backpropation step"""
                while node!= None:
                    #print "node.parentNode:",node.parentNode
                    if re <= 0:
                        node.Update(re, add_vis_count = 0)
                    else:
                        node.Update(re, add_vis_count = 1)
                    node.delete_virtual_loss()
                    node = node.parentNode
                lock.release()


    result_queue.put([all_compounds,wave_compounds,depth,wave,maxnum,uv_intensity_list,deen_list,gap_list,reward_list,index_list,mol_weight_list,logP_list,SA_score_list,depth_list])

def charge_check(mol):
    print 'charge_checking'
    standard_valence_list = [0, 1, 0, 1, 2, 3, 4, 3, 2, 1, 0, 1, 2, 3, 4, 3, 2, 1, 0, 1, 2]
    check = True
    for atom in mol.GetAtoms():
        if standard_valence_list[atom.GetAtomicNum()] != atom.GetExplicitValence():
            check = False
            break
    return check

def gaussion_workers(chem_model,val,gau_parallel,charge_check):
    while True:
        simulation_time=time.time()
        task = comm.recv(source=0, tag=MPI.ANY_TAG, status=status)
        tag = status.Get_tag()
        if tag==START:
            state=task[0]
            m=task[1]
            ind=task[2]
            all_posible=chem_kn_simulation(chem_model,state,val,m)
            generate_smile=predict_smile(all_posible,val)
            new_compound=make_input_smile(generate_smile)
            score=[]
            kao=[]
            intensity = -1000000
            deen = 1000000
            gap = 1000000
            mol_weight = 0
            SA_score = 10
            logP = 0
            dp = len(state)
            intensity_list = []
            wl_list = []
            standard_valence_list = [0, 1, 0, 1, 2, 3, 4, 3, 2, 1, 0, 1, 2, 3, 4, 3, 2, 1, 0, 1, 2]
            s1_wavenum=-1000
            s1_strength=0
            s1_wavelength = None
            s1_wl_list = []
            s1_str_list = []

            try:
                m = Chem.MolFromSmiles(str(new_compound[0]))
                mol_weight = Descriptors.MolWt(m)
                logP = Crippen.MolLogP(m)
                SA_score = sascorer.calculateScore(m)
                #print 'prev add Hs'
                m_H = Chem.AddHs(m)
                #print Chem.MolToSmiles(m_H)
                
                standard_valence_list = [0, 1, 0, 1, 2, 3, 4, 3, 2, 1, 0, 1, 2, 3, 4, 3, 2, 1, 0, 1, 2]
                ccheck = True
                if charge_check:
                    for atom in m_H.GetAtoms():
                        #print standard_valence_list[atom.GetAtomicNum()], atom.GetExplicitValence()
                        if standard_valence_list[atom.GetAtomicNum()] != atom.GetExplicitValence():
                            ccheck = False
                            break

                if not ccheck:
                    #print 'charge checked: not'
                    m = None
                #print 'charge checked'
                #if charge_check(m_H, standard_valence_list):
                #    print('charge_check: OK')
                #else:
                #    m = None

            except:
                m=None
            #if m!=None and len(task[i])<=81:


            if m!=None:
                try:
                    stable=tansfersdf(str(new_compound[0]),ind)
                except:
                    stable = -1
                if stable==1.0:
                    cd_path = os.getcwd()
                    try:
			SDFinput = 'CheckMolopt'+str(ind)+'.sdf'
                        #wavelength=GauTDDFT_ForDFT('B3LYP', '3-21G*', 1, 'CheckMolopt'+str(ind)+'.sdf')
                        calc_sdf = GaussianDFTRun('B3LYP', '3-21G*', gau_parallel, 'OPT fluor energy deen uv homolumo', SDFinput, 0)
                        outdic = calc_sdf.run_gaussian()
                        wavelength = outdic['uv'][0]
                        s1_wavelength = outdic['S1 Wavelength and Oscillator strengths'][0]

                        if os.path.isfile('CheckMol'+str(ind)+'.sdf'):
                            shutil.move('CheckMol'+str(ind)+'.sdf', 'dft_result')
                        if os.path.isfile('CheckMolopt'+str(ind)+'.sdf'):
                            shutil.move('CheckMolopt'+str(ind)+'.sdf', 'dft_result')
                        #if os.path.isfile('CheckMolopt'+str(ind)+'.com'):
                        #    shutil.move('CheckMolopt'+str(ind)+'.com', 'dft_result')
                        #if os.path.isfile('CheckMolopt'+str(ind)+'.log'):
                        #    shutil.move('CheckMolopt'+str(ind)+'.log', 'dft_result')
                        #if os.path.isfile('CheckMolopt'+str(ind)+'.chk'):
                        #    os.remove('CheckMolopt'+str(ind)+'.chk')
                    except:
                        os.chdir(cd_path)
                        wavelength=None
                        if os.path.isfile('CheckMolopt'+str(ind)+'.sdf'):
                            os.remove('CheckMolopt'+str(ind)+'.sdf')
                        if os.path.isfile('CheckMol'+str(ind)+'.sdf'):
                            os.remove('CheckMol'+str(ind)+'.sdf')
                else:
                    wavelength=None
                    if os.path.isfile('CheckMolopt'+str(ind)+'.sdf'):
                        os.remove('CheckMolopt'+str(ind)+'.sdf')
                    if os.path.isfile('CheckMol'+str(ind)+'.sdf'):
                        os.remove('CheckMol'+str(ind)+'.sdf')

                if wavelength!=None and wavelength!=[]:
                    wavenum=wavelength[0]
                    intensity=outdic['uv'][1][0]
                    deen=outdic['deen']
                    gap=outdic['gap']
                    wl_list = outdic['uv'][0]
                    intensity_list = outdic['uv'][1]
                else:
                    wavenum=-1000
                
                if s1_wavelength!=None and s1_wavelength != []:
                    s1_wavenum=s1_wavelength[0]
                    s1_strength=outdic['S1 Wavelength and Oscillator strengths'][1][0]
                    s1_wl_list=outdic['S1 Wavelength and Oscillator strengths'][0]
                    s1_str_list=outdic['S1 Wavelength and Oscillator strengths'][1]
                else:
                    s1_wavenum=-1000
            else:
                wavenum=-1000
            score.append(wavenum)
            score.append(new_compound[0])
            score.append(rank)
            score.append(intensity)
            score.append(deen)
            score.append(gap)
            score.append(wl_list)
            score.append(intensity_list)
            score.append(ind)
            score.append(mol_weight)
            score.append(logP)
            score.append(SA_score)
            score.append(dp)
            score.append(s1_wavenum) #flour 
            score.append(s1_strength)
            score.append(s1_wl_list)
            score.append(s1_str_list)

            comm.send(score, dest=0, tag=DONE)
            simulation_fi_time=time.time()-simulation_time
            print("simulation_fi_time:",simulation_fi_time)
        if tag==EXIT:
            MPI.Abort(MPI.COMM_WORLD)


    comm.send([-1000,'',0,0,0,0,[],[],ind,0,0,0,0, -1000, 0, [], []], dest=0, tag=EXIT)


if __name__ == "__main__":
    comm=MPI.COMM_WORLD
    size=comm.size
    rank=comm.rank
    status=MPI.Status()
    READY, START, DONE, EXIT = 0, 1, 2, 3

    #smile_old=zinc_data_with_bracket_original()
    #val,smile=zinc_processed_with_bracket(smile_old)
    #val=['\n', '&', 'C', '[C@@H]', '(', 'N', ')', 'O', '=', '1', '/', 'c', 'n', '[nH]', '[C@H]', '2', '[NH]', '[C]', '[CH]', '[N]', '[C@@]', '[C@]', 'o', '[O]', '3', '#', '[O-]', '[n+]', '[N+]', '[CH2]', '[n]'] #original
    #atoms zinc 250,000
    
    #val=['\n', '&', 'C', '(', ')', 'c', '1', '2', 'o', '=', 'O', 'N', '3', 'F', '[C@@H]', 'n', '-', '#', 'S', 'Cl', '[O-]', '[C@H]', '[NH+]', '[C@]', 's', 'Br', '/', '[nH]', '[NH3+]', '4', '[NH2+]', '[C@@]', '[N+]', '[nH+]', '\\', '[S@]', '5', '[N-]', '[n+]', '[S@@]', '[S-]', '6', '7', 'I', '[n-]', 'P', '[OH+]', '[NH-]', '[P@@H]', '[P@@]', '[PH2]', '[P@]', '[P+]', '[S+]', '[o+]', '[CH2-]', '[CH-]', '[SH+]', '[O+]', '[s+]', '[PH+]', '[PH]', '8', '[S@@+]']
    val=['\n', '&', 'C', '(', ')', 'c', '1', '2', 'o', '=', 'O', 'N', '3', 'F', '[C@@H]', 'n', '-', '#', '/', '[nH]', 'Br', '[C@H]', 'Cl', '[C@]', '[C@@]', '\\', '4', '5', '6', '7', 'I']

    chem_model=loaded_model()
    graph = tf.get_default_graph()
    chemical_state = chemical()

    ts_strategy = 'puct' #'uct', 'puct' 
    search_parameter = 0.25 #If ts_strategy=='uct', 0 < search_parameter < 1. If ts_strategy=='puct', default value is 5 (AlphaGo). 
    num_simulations = 29 # core - 1, max: 2560 (skylake)
    gau_parallel = 1
    num_rollout = 3
    simulation_time = 3600*120 # 3600*24 # max: 168h
    alpha = 1 # alph*mean + (1 - alpha)*max + bais
    objective = 'WLFL' # 'WL_IT', 'HL', 'WL'
    charge_check = True # True or False
    SA_score_check = True # True or False
    output_file = 'csvresult_FP_'+ts_strategy+'_C'+str(search_parameter)+'_alpha'+str(alpha)+'_obj'+objective+'_para'+str(num_simulations)+'_time'+str(simulation_time/3600)+'h_rollout'+str(num_rollout)+'_CC'+str(charge_check)+'_SA'+str(SA_score_check)+'_0218.csv'

    thread_pool=[]
    lock=Lock()
    gau_file_index=0

    """initialization of the chemical trees and grammar trees"""
    root=['&']
    rootnode = Node(position= root)
    maxnum=0
    ind_mol=0
    reward_dis=[]
    all_compounds=[]
    wave_compounds=[]
    #result_test=[]
    wave=[]
    deen_list = []
    gap_list = []
    uv_intensity_list = []
    wl_list_list = []
    intensity_list_list = []
    reward_list = []
    index_list = []
    mol_weight_list = []
    logP_list = []
    SA_score_list = []
    depth_list = []
    s1_wavelength_list = []
    s1_strength_list = []
    s1_wavelength_list_list = []
    s1_strength_list_list = []
    depth=[]
    result=[]
    result_queue=Queue()
    free_core_id=range(1,num_simulations+1)
    use_core_id = []
    #free_core_id=[40,80,120,160,200,240,280,320,360,400,440,480,520,560,600,640,680,720,760,800,840,880,920,960,1000]

    if rank==0:
        for thread_id in range(num_simulations):
            thread_best = Thread(target=ChemTS_run,args=(rootnode,result_queue,lock,chem_model,ts_strategy,search_parameter,num_simulations, gau_parallel,simulation_time,output_file,alpha,objective,num_rollout,charge_check,SA_score_check))
            thread_pool.append(thread_best)

        for i in range(num_simulations):
            thread_pool[i].start()

        for i in range(num_simulations):
            thread_pool[i].join()
        for i in range(num_simulations):
            result.append(result_queue.get())
        #print(result[0][1])
        #print(result[0][2])
        #print(result[0][3])
        #print(len(result[0][0]))
        #print(len(result[0][1]))
        #print(max(result[0][3]))
        #print(len(result[0][3]))
        #print(result[0][4])
        comm.Abort()
        for i in range(len(free_core_id)):
            comm.send(None, dest=i+1, tag=EXIT)
    #elif rank%40==0 and rank!=0:
    else:
	    gaussion_workers(chem_model,val,gau_parallel, charge_check)
    #else:
       # while True:
           # time.sleep(30)
           
