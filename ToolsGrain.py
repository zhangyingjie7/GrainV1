#-------------------------------------------------------------------------------
# Name:        Tools
# Purpose:     A framework for fast correlation analysis
#
# Author:      Zhang Yingjie
#
# Created:     12-07-2019
# Copyright:   (c)  2019
#
# Usage: from Tools import *
#
# Class Basic_Tools: Some tools to process Term and variables
#
# Class Correlation_of_Term_F: Some tools to compute correlation of Term or F 
#
# Class LinearApproximate_of_H: To compute linear approximate representation of H function
#
#-------------------------------------------------------------------------------
from sympy import *
from sympy.logic import simplify_logic  # CNF or DNF

from itertools import combinations

import networkx as nx
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings("ignore")  # ignore warning

import numpy as np
import math
import time


global flag # if ['1'] in Term, flag=1; else flag=0

class Basic_Tools:

    @staticmethod
    def wordToBinaryString(word, size):
        """
        Example
        ----------
        Basic_Tools.wordToBinaryString(0xF1, 8) >>>'11110001'
        """
        return bin(word)[2:].zfill(size)
    
    def binaryStringToWord(self,string):
        word=0
        for i in range(0,len(string)):
            word=word^(int(string[i])<<(len(string)-1-i))
        return word
        
    
    def reduce_Term(self,Term):
        """
        Firstly, sort Term; then
        if term appears even times in Term, remove term;
        if term appears odd times in Term, reserve term.
        
        Example
        ----------
        Term=[['b_102', 'b_104'],['b_102', 'b_104'],['b_102', 'b_104'],['b_101', 'b_104'],['b_101', 'b_104']]
        Term=reduce_Term(Term)
        # >>>Term
        [['b_102', 'b_104']]
        """
        global flag
        if ['1'] in Term:
            flag=1
            Term.remove(['1'])
        else:
            flag=0
        # print('flag',flag)    
        Term=sorted(Term)  # Sort Term in an order, for the aim of removing correctly.
        i=len(Term)-1
        while i>=1:
            if Term[i]==Term[i-1]:
                Term.remove(Term[i])
                Term.remove(Term[i-1])
                i=i-2
            else:
                i=i-1
        return Term
    
    def genVars_From_Term(self,Term): 
        """
        Generate variable list from term
        
        Example
        ----------  
        T=Basic_Tools()
        Term=[['n_27'], ['n_27', 'n_33'], ['n_34']]
        V=T.genVars_From_Term(Term)
        for v in V:print(v,end=' ')
        # >>> n_27 n_33 n_34
        """
        Vars=set()
        for i in range(0,len(Term)):
            Vars=Vars|set(Term[i])
        return sorted(list(Vars))
    
    def transTerm_to_ANF(self,Term): 
        """
        Transform term to its ANF
        
        Example
        ---------- 
        T=Basic_Tools()
        Term=[['n_27'], ['n_27', 'n_33'], ['n_34']]
        f=T.transTerm_to_ANF(Term)
        print(f)
        #>>>
        Xor(n_27, n_34, n_27 & n_33)
        """
        temp_str=Term[0][0]
        for j in range(1,len(Term[0])):
            temp_str=temp_str+ ' & ' + Term[0][j] 
        f=temp_str
        for i in range(1,len(Term)):
            temp_str=Term[i][0]
            for j in range(1,len(Term[i])):
                temp_str=temp_str+ ' & ' + Term[i][j]
            f=f + ' ^ ' + temp_str
        """
        f=f+' ^ '+str(1)
        print(f)
        """
        f=sympify(f,convert_xor=False)
        return f
    
    def transTerm_to_graphFormat(self,Term):
        """
        Transform term to Latex format
        
        Example
        ----------
        T=Basic_Tools()
        Term=[['b_102', 'b_104']]
        Term=T.transTerm_to_graphFormat(Term)
        # >>>Term
        ['$b_{102}$', '$b_{104}$']
        
        """
        for i in range(0,len(Term)):
            for j in range(0,len(Term[i])):
                Term[i][j]='$'+Term[i][j][0:2]+'{'+Term[i][j][2:]+'}$'
        return Term
    
    def transGraph_to_term(self,Graph_Term):
        """
        Transform Latex format to term, then term can be transformed into ANF format
        
        Example
        ----------
        T=Basic_Tools()
        Term=['$b_{102}$', '$b_{104}$']
        Term=T.transGraph_to_term(Term)
        # >>>Term
        [['b_102', 'b_104']]
        """
        for i in range(0,len(Graph_Term)):
            for j in range(0,len(Graph_Term[i])):
                Graph_Term[i][j]=Graph_Term[i][j].replace('$','')
                Graph_Term[i][j]=Graph_Term[i][j].replace('{','')
                Graph_Term[i][j]=Graph_Term[i][j].replace('}','')
        return Graph_Term
    
    @staticmethod
    def H_corr_dict(H_corr_value):
        """
        Parameters
        ----------
        H_corr_value: list of H_corr
        
        Returns
        ----------
        dict_list: Nonzero dictionary of H_corr_value
        H_index: Index list of nonzero corr in H_corr_value
        
        Example
        ----------
        H_corr_value=[[32, 0, 32, 0, 32, 0, -32, 0, 32, 0, 32, 0, 32, 0, -32, 0, 32, 0, 32, 0, 32, 0, -32, 0, 32, 0, 32, 0, 32, 0, -32, 0, 32, 0, 32, 0, 32, 0, -32, 0, 32, 0, 32, 0, 32, 0, -32, 0, -32, 0, -32, 0, -32, 0, 32, 0, -32, 0, -32, 0, -32, 0, 32, 0, 32, 0, 32, 0, 32, 0, -32, 0, 0, 32, 0, 32, 0, 32, 0, -32, 32, 0, 32, 0, 32, 0, -32, 0, 0, 32, 0, 32, 0, 32, 0, -32, 32, 0, 32, 0, 32, 0, -32, 0, 0, 32, 0, 32, 0, 32, 0, -32, -32, 0, -32, 0, -32, 0, 32, 0, 0, -32, 0, -32, 0, -32, 0, 32],
            [32, 0, 32, 0, 32, 0, -32, 0, -32, 0, -32, 0, -32, 0, 32, 0, 32, 0, 32, 0, 32, 0, -32, 0, -32, 0, -32, 0, -32, 0, 32, 0, 32, 0, 32, 0, 32, 0, -32, 0, -32, 0, -32, 0, -32, 0, 32, 0, -32, 0, -32, 0, -32, 0, 32, 0, 32, 0, 32, 0, 32, 0, -32, 0, 32, 0, 32, 0, 32, 0, -32, 0, 0, -32, 0, -32, 0, -32, 0, 32, 32, 0, 32, 0, 32, 0, -32, 0, 0, -32, 0, -32, 0, -32, 0, 32, 32, 0, 32, 0, 32, 0, -32, 0, 0, -32, 0, -32, 0, -32, 0, 32, -32, 0, -32, 0, -32, 0, 32, 0, 0, 32, 0, 32, 0, 32, 0, -32],
            [32, 0, 32, 0, 32, 0, -32, 0, 32, 0, 32, 0, 32, 0, -32, 0, 32, 0, 32, 0, 32, 0, -32, 0, 32, 0, 32, 0, 32, 0, -32, 0, 32, 0, 32, 0, 32, 0, -32, 0, 32, 0, 32, 0, 32, 0, -32, 0, -32, 0, -32, 0, -32, 0, 32, 0, -32, 0, -32, 0, -32, 0, 32, 0, -32, 0, -32, 0, -32, 0, 32, 0, 0, -32, 0, -32, 0, -32, 0, 32, -32, 0, -32, 0, -32, 0, 32, 0, 0, -32, 0, -32, 0, -32, 0, 32, -32, 0, -32, 0, -32, 0, 32, 0, 0, -32, 0, -32, 0, -32, 0, 32, 32, 0, 32, 0, 32, 0, -32, 0, 0, 32, 0, 32, 0, 32, 0, -32],
            [32, 0, 32, 0, 32, 0, -32, 0, -32, 0, -32, 0, -32, 0, 32, 0, 32, 0, 32, 0, 32, 0, -32, 0, -32, 0, -32, 0, -32, 0, 32, 0, 32, 0, 32, 0, 32, 0, -32, 0, -32, 0, -32, 0, -32, 0, 32, 0, -32, 0, -32, 0, -32, 0, 32, 0, 32, 0, 32, 0, 32, 0, -32, 0, -32, 0, -32, 0, -32, 0, 32, 0, 0, 32, 0, 32, 0, 32, 0, -32, -32, 0, -32, 0, -32, 0, 32, 0, 0, 32, 0, 32, 0, 32, 0, -32, -32, 0, -32, 0, -32, 0, 32, 0, 0, 32, 0, 32, 0, 32, 0, -32, 32, 0, 32, 0, 32, 0, -32, 0, 0, -32, 0, -32, 0, -32, 0, 32]]  
        H_corr_dict,H_index=Basic_Tools.H_corr_dict(H_corr_value)
        
        for v in H_corr_dict:
            print(v)
        
        print()
        for v in H_index:
            print(v)
        # >>>
        {0: 32, 2: 32, 4: 32, 6: -32, 8: 32, 10: 32, 12: 32, 14: -32, 16: 32, 18: 32, 20: 32, 22: -32, 24: 32, 26: 32, 28: 32, 30: -32, 32: 32, 34: 32, 36: 32, 38: -32, 40: 32, 42: 32, 44: 32, 46: -32, 48: -32, 50: -32, 52: -32, 54: 32, 56: -32, 58: -32, 60: -32, 62: 32, 64: 32, 66: 32, 68: 32, 70: -32, 73: 32, 75: 32, 77: 32, 79: -32, 80: 32, 82: 32, 84: 32, 86: -32, 89: 32, 91: 32, 93: 32, 95: -32, 96: 32, 98: 32, 100: 32, 102: -32, 105: 32, 107: 32, 109: 32, 111: -32, 112: -32, 114: -32, 116: -32, 118: 32, 121: -32, 123: -32, 125: -32, 127: 32}
        {0: 32, 2: 32, 4: 32, 6: -32, 8: -32, 10: -32, 12: -32, 14: 32, 16: 32, 18: 32, 20: 32, 22: -32, 24: -32, 26: -32, 28: -32, 30: 32, 32: 32, 34: 32, 36: 32, 38: -32, 40: -32, 42: -32, 44: -32, 46: 32, 48: -32, 50: -32, 52: -32, 54: 32, 56: 32, 58: 32, 60: 32, 62: -32, 64: 32, 66: 32, 68: 32, 70: -32, 73: -32, 75: -32, 77: -32, 79: 32, 80: 32, 82: 32, 84: 32, 86: -32, 89: -32, 91: -32, 93: -32, 95: 32, 96: 32, 98: 32, 100: 32, 102: -32, 105: -32, 107: -32, 109: -32, 111: 32, 112: -32, 114: -32, 116: -32, 118: 32, 121: 32, 123: 32, 125: 32, 127: -32}
        {0: 32, 2: 32, 4: 32, 6: -32, 8: 32, 10: 32, 12: 32, 14: -32, 16: 32, 18: 32, 20: 32, 22: -32, 24: 32, 26: 32, 28: 32, 30: -32, 32: 32, 34: 32, 36: 32, 38: -32, 40: 32, 42: 32, 44: 32, 46: -32, 48: -32, 50: -32, 52: -32, 54: 32, 56: -32, 58: -32, 60: -32, 62: 32, 64: -32, 66: -32, 68: -32, 70: 32, 73: -32, 75: -32, 77: -32, 79: 32, 80: -32, 82: -32, 84: -32, 86: 32, 89: -32, 91: -32, 93: -32, 95: 32, 96: -32, 98: -32, 100: -32, 102: 32, 105: -32, 107: -32, 109: -32, 111: 32, 112: 32, 114: 32, 116: 32, 118: -32, 121: 32, 123: 32, 125: 32, 127: -32}
        {0: 32, 2: 32, 4: 32, 6: -32, 8: -32, 10: -32, 12: -32, 14: 32, 16: 32, 18: 32, 20: 32, 22: -32, 24: -32, 26: -32, 28: -32, 30: 32, 32: 32, 34: 32, 36: 32, 38: -32, 40: -32, 42: -32, 44: -32, 46: 32, 48: -32, 50: -32, 52: -32, 54: 32, 56: 32, 58: 32, 60: 32, 62: -32, 64: -32, 66: -32, 68: -32, 70: 32, 73: 32, 75: 32, 77: 32, 79: -32, 80: -32, 82: -32, 84: -32, 86: 32, 89: 32, 91: 32, 93: 32, 95: -32, 96: -32, 98: -32, 100: -32, 102: 32, 105: 32, 107: 32, 109: 32, 111: -32, 112: 32, 114: 32, 116: 32, 118: -32, 121: -32, 123: -32, 125: -32, 127: 32}
        
        [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62, 64, 66, 68, 70, 73, 75, 77, 79, 80, 82, 84, 86, 89, 91, 93, 95, 96, 98, 100, 102, 105, 107, 109, 111, 112, 114, 116, 118, 121, 123, 125, 127]
        [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62, 64, 66, 68, 70, 73, 75, 77, 79, 80, 82, 84, 86, 89, 91, 93, 95, 96, 98, 100, 102, 105, 107, 109, 111, 112, 114, 116, 118, 121, 123, 125, 127]
        [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62, 64, 66, 68, 70, 73, 75, 77, 79, 80, 82, 84, 86, 89, 91, 93, 95, 96, 98, 100, 102, 105, 107, 109, 111, 112, 114, 116, 118, 121, 123, 125, 127]
        [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62, 64, 66, 68, 70, 73, 75, 77, 79, 80, 82, 84, 86, 89, 91, 93, 95, 96, 98, 100, 102, 105, 107, 109, 111, 112, 114, 116, 118, 121, 123, 125, 127]

        """

        dict_list=[{} for i in range(0,len(H_corr_value))]
        for i in range(0,len(H_corr_value)):
            for j in range(0,len(H_corr_value[0])):
                if H_corr_value[i][j]!=0:
                    dict_list[i].update({j:H_corr_value[i][j]})
        H_index=[[] for i in range(0,len(dict_list))]
        for i in range(0,len(H_index)):
            for k in dict_list[i].keys():
                H_index[i].append(k)
        return dict_list,H_index
    
    @staticmethod
    def gen_SpanW(alpha,start_point):
        W=[]
        spanLen=len(alpha)
        for i in range(0,2**spanLen):
            temp=[0 for i in range(0,len(alpha[0]))]
            index=list(map(int, Basic_Tools.wordToBinaryString(i, spanLen)))
            for j in range(0,len(alpha)):
                if index[j]==1:
                    temp=list(map(lambda x: x[0]^x[1],zip(alpha[j]*index[j],temp)))
            temp=list(map(lambda x: x[0]^x[1],zip(start_point,temp)))
            W=W+[temp]
        return W

    
class LinearApproximate_of_H(Basic_Tools):
    
    def __init__(self,Term,b_name,file):
        """
        Parameters
        ----------
        Term: term list of ANF
        b_name: variables of NFSR
        file: file to save the approximation
        """
        self.Term=Term
        self.b_name=b_name
        self.filename=file
    
    def genValue_of_productSum(self, a, x):
        """
        Parameters
        ----------
        a: mask
        x: values
        
        Returns
        ----------
        product: (a,x)
        """
        a_x = 0
        for i in range(0,len(a)):
            a_x = a_x ^ (int (a[i]) * int(x[i]))
        return a_x
    
    def get_sub_approximate(self,k,b_mask):
        """
        Parameters
        ----------
        k: mask of other variables
        b_mask: mask of b variables
        
        Notice
        ----------
        1. To keep variables in order, we need to write 's_08' rather than 's_8'
        2. To make sure b_name is in front of s_name, change some details if it isn't.
        """
        f=self.transTerm_to_ANF(self.Term)
        x=self.genVars_From_Term(self.Term)
        assert k<=2**(len(x)-len(self.b_name)), "Wrong linear mask!!!"
        mask=list(self.wordToBinaryString(k,(len(x)-len(self.b_name))))
        count0 = 0
        count1 = 0
        for i in range(0,2**len(x)):
            f_temp=f
            x_value=list(self.wordToBinaryString(i, len(x)))
            for j in range(0,len(x_value)):
                f_temp=f_temp.subs([(symbols(x[j]),int(x_value[j]))])
            # f_temp's value: {0,1,True,False}
            if f_temp==True:
                f_temp=1
            if f_temp==False:
                f_temp=0
            if self.genValue_of_productSum((b_mask+mask),x_value) == f_temp:
            # if b_name is in the back, then mask+b_mask
                count0 = count0 + 1
            else:
                count1 = count1 + 1
        return (count0 - count1)
    
    def get_approximate(self):
        """
        Returns
        corr_list: Mask table of Term divided by b_value
        
        """
        b_name=self.b_name
        x=self.genVars_From_Term(self.Term)
        corr_list=[[] for i in range(0,2**len(b_name))]
        for i in range(0,2**len(b_name)):
            b_value=list(self.wordToBinaryString(i, len(b_name)))
            print(b_value,file=self.filename)
            for k in range(0,2**(len(x)-len(self.b_name))):
                corrTemp=self.get_sub_approximate(k, b_value)
                corr_list[i]=corr_list[i]+[corrTemp]
            print(corr_list[i])
            print(corr_list[i],file=self.filename)
            print(file=self.filename)
        return corr_list
        
        
        
class Graph(Basic_Tools):
    """
    Functions
    ----------
    1. genGraph(Term): Generate graph of a Term
    2. showGraph(G): Show picture of G
    3. kCutSets_of_graph(G,k): Find optimal k node cut-sets of G, see details for meaning of "optimal"
    4. find_validCutSet_of_graph(G): Find a valid and optimal cut-sets of G
    5. cutSets_of_Graph(G): Find a cut-sets of G, the number of vertices of each Sub-graphs <= Threshold
    """

    def __init__(self,Threshold,k):
        """
        Parameters
        ----------
        Threshold: maximum number of vertices in Sub-graphs
        """
        
        self.Threshold=Threshold
        self.Init_k=k

    
    def genGraph(self,Term):
        """
        Generate graph of a Term
        Graph nodes: variables of Term
        Graph edges: generate by every single term 
        
        Example
        ----------
        Term=[['n_27'], ['n_27', 'n_33'], ['n_34']]
        Nodes: 'n_27','n_33','n_34'
        Edges: (n_27,n_33)
        
        Notice
        ----------
        To show graph clearly, we need to transform Term into Latex format:
        Term=Basic_Tools.transTerm_to_graphFormat(Term)
        """
        Vars=self.genVars_From_Term(Term) 
        edges=[]
        for v in Term:
            t=len(v)
            if t==1:
                edges=edges+[]
            else:
                temp=[v[i] for i in range(0,t)]
                edges=edges+list(combinations(temp,2))
        
        G=nx.Graph()
        G.add_nodes_from(Vars)
        G.add_edges_from(edges)
        
        '''
        # show figure
        plt.figure(1)
        nx.draw(G, with_labels = True,font_size=8,node_size=200,width=1)
        plt.draw()
        # plt.savefig('Figure_1.png',dpi=100)  # Not clear, save manually.
        plt.show() # close the figure or remove this statement, program will run the following statement
        '''
        return G
    
    def showGraph(self,G):
        """
        Parameters
        ----------
        G: Graph
        
        Note
        ----------
        1. Save picture manually to ensure its quality.
        2. Close the window after it is saved manually.
        """
        plt.figure()
        nx.draw(G, with_labels = True,font_size=8,node_size=200,width=1)
        plt.draw()
        plt.show()
    
    def kCutSets_of_graph(self,G,k):
        """
        Parameters
        ----------
        G: Graph
        k: To find optimal k node cut-sets of G
        
        Returns
        ----------
        G_cutsets: Optimal k cut-sets of G,
        which are to removed to make most sub-graphs with nodes less than Threshold 
        """
        Threshold=self.Threshold
        G_subs=list(nx.connected_components(G)) # all connected components of G
        assert min([len(list(c)) for c in G_subs]) != 1, 'Graph with single node, zero correlation'
        G_cutsets=[] # cut points of G 
        for c in G_subs:
            if len(list(c))>Threshold:
                G1=G.subgraph(c).copy()               
                cutsets = list(nx.all_node_cuts(G1, k))  # all cut_size node cuts  
                if cutsets==[]:
                    return []
                arr_var=10000000
                for i in range(0,len(cutsets)):
                    G2=G1.copy()
                    G2.remove_nodes_from(cutsets[i])
                    G2_subs=list(nx.connected_components(G2))
                    
                    
                    subs_len=[len(v) for v in G2_subs]
                    if np.var(subs_len)< arr_var: # To find cut with minimum variance
                    # var = mean(abs(x - x.mean())**2)
                        cut=cutsets[i]
                        arr_var=np.var(subs_len)
                    '''
                    subs_len=[abs(len(v)-self.Threshold) for v in G2_subs]
                    if sum(subs_len)< arr_var: # To find cut which is the most closest to self.Threshold
                        cut=cutsets[i]
                        arr_var=sum(subs_len) 
                    '''                     
                G_cutsets=G_cutsets+list(cut)   
        return G_cutsets
    
    def find_validCutSet_of_graph(self,G):
        
        """
        Parameters
        ----------
        G: Graph
        
        Returns
        ----------
        cut_sets: A valid optimal cut-sets of G
        - Minimal number of nodes
        - Minimal variance of sub-graphs' length
        """
        k=0
        cut_sets=[]
        while cut_sets==[]: # To make sure "cut_sets" isn't empty
            k=k+1
            cut_sets=self.kCutSets_of_graph(G, k) # Optimal k cut-sets of G
        return cut_sets
    
    def cutSets_of_Graph(self,G):
        """
        Parameters
        ----------
        G: Graph
        
        Note
        ----------
        Show graph, then save and close.
        
        Returns
        ----------
        cut_sets: A cut-sets of G, the number of vertices of each Sub-graphs <= Threshold
        
        """
        Threshold=self.Threshold
        
        # 1st cut sets of G
        if self.Init_k==1:
            G_cutsets=self.find_validCutSet_of_graph(G)
        else:
            G_cutsets=self.kCutSets_of_graph(G, self.Init_k)
        
        G1=G.copy() # Copy of G
        G1.remove_nodes_from(G_cutsets)
        subGraph = list(nx.connected_components(G1))  
        sub_len=[]
        for v in subGraph:
            sub_len=sub_len+[len(v)]
          
        while max(sub_len) > Threshold:
            for v in subGraph: # type(v)=<class 'set'>
                v=list(v)
                if len(v)>Threshold:
                    G2=G1.subgraph(v).copy()
                    G_cutsets=G_cutsets + self.find_validCutSet_of_graph(G2)
            G1=G.copy() # A new copy of G
            G1.remove_nodes_from(G_cutsets)
            subGraph = list(nx.connected_components(G1))  
            sub_len=[]
            for v in subGraph:
                sub_len=sub_len+[len(v)]
        '''
        G1=G.copy() # Copy of G
        G1.remove_nodes_from(G_cutsets)
        self.showGraph(G1)
        '''
        return G_cutsets
    
    def divide_Term(self,Term):
        G=self.genGraph(Term)
        sub = list(sorted(nx.connected_components(G),key=len,reverse=False)) 
        
        # Connected components of G
        sub_len=[]
        sub_vars=[]
        for v in sub:
            v=list(v)
            sub_len=sub_len+[len(v)]
            sub_vars=sub_vars+[set(v)]
        G_cutsets=[]
        
        # The maximum size of G_connected_components is larger than self.Threshold
        if max(sub_len)>self.Threshold:
            G_cutsets=self.cutSets_of_Graph(G)
            G.remove_nodes_from(G_cutsets) 
            sub = list(sorted(nx.connected_components(G),key=len,reverse=False))
            sub_len=[]
            sub_vars=[]
            for v in sub:
                v=list(v)
                sub_len=sub_len+[len(v)]
                sub_vars=sub_vars+[set(v)]
            
        divided_Term=[[] for i in range(0,len(sub_vars))]
        for i in range(0,len(sub_vars)):
            for v in range(0,len(Term)):
                if set(Term[v]) & sub_vars[i] != set():
                    divided_Term[i]=divided_Term[i]+[Term[v]]
        
        Term_merge=[]
        for v in divided_Term:
            Term_merge=Term_merge+v
            
        for v in Term:
            if v not in Term_merge:
                divided_Term[0]=divided_Term[0]+[v]  # add cut-points to divided_Term[0]
    
        return G_cutsets,divided_Term

  
class Correlation_of_Term_F(Graph):
    """
    To compute correlation of a term 
    Term=[['n_27'], ['n_27', 'n_33'], ['n_34']] # means: n27 ^ n27*n33 ^ n34
    
    """
    def correlation_of_smallTerm(self,Term):
        """
        Directly compute correlation of Term with <= Threshold variables
        
        Example
        ---------- 
        T=Correlation_of_Term_F(10)
        Term=[['n_27'], ['n_27', 'n_33']]
        F=T.Correlation_of_Term(Term)
        print(F)
        # >>>
        0.5
        """
        
        x=self.genVars_From_Term(Term)
        f=self.transTerm_to_ANF(Term)
        count0=0
        count1=0
        for k in range(0,2**len(x)):
            x_values=list(self.wordToBinaryString(k, len(x)))
            #x_values=list(bin(k)[2:].zfill(len(x)))
            f_temp=f
            for x_index in range(0,len(x)):
                f_temp=f_temp.subs([(symbols(x[x_index]),int(x_values[x_index]))])
            if f_temp==Integer(1):
                count1=count1+1 
            elif f_temp==True:
                count1=count1+1
            elif f_temp==Integer(0):
                count0=count0+1
            elif f_temp==False:
                count0=count0+1
        return (count0-count1)/(2**len(x))
    
    def correlation_of_smallF(self,f,x):
        """
        Compute correlation of Term
        
        Parameters
        ----------
        f: ANF of term
        x: list of variables in f
        
        Example
        ----------
        T=Correlation_of_Term_F(10)
        Term=[['n_27'], ['n_27', 'n_33']]
        f=T.transTerm_to_ANF(Term)
        x=T.genVars_From_Term(Term)
        print(T.correlation_of_smallF(f,x))
        # >>>
        0.5
        """
        count0=0
        count1=0
        for k in range(0,2**len(x)):
            #x_values=list(bin(k)[2:].zfill(len(x)))
            x_values=list(self.wordToBinaryString(k, len(x)))
            f_temp=f
            for x_index in range(0,len(x)):
                f_temp=f_temp.subs([(symbols(x[x_index]),int(x_values[x_index]))])
            if f_temp==Integer(1):
                count1=count1+1 
            elif f_temp==True:
                count1=count1+1
            elif f_temp==Integer(0):
                count0=count0+1
            elif f_temp==False:
                count0=count0+1
        return (count0-count1)/(2**len(x))
    
    def correlation_of_Term_with_cutpoints(self,Term,Cut_points,myfile):
        Vars=self.genVars_From_Term(Term)
        f=self.transTerm_to_ANF(Term)
        corr_list=[0 for i in range(0,2**len(Cut_points))]
        if set(Vars) & set(Cut_points)!=set():
            for cut_len in range(0,2**len(Cut_points)):
                #cut_values=list(bin(cut_len)[2:].zfill(len(Cut_points)))
                cut_values=list(self.wordToBinaryString(cut_len, len(Cut_points)))
                f_tem=f
                for k in range(0,len(Cut_points)):
                    f_tem=f_tem.subs([(symbols(Cut_points[k]),int(cut_values[k]))])
                corr_list[cut_len]=self.correlation_of_smallF(f_tem, list(set(Vars)-set(Cut_points))) 
            print(Term,corr_list,file=myfile)
            return corr_list
        else:
            corr=self.correlation_of_smallF(f, Vars)
            print(Term,[corr],file=myfile)
            return [corr]
        
    def correlation_of_Term_with_cutpointsNonFile(self,Term,Cut_points):
        Vars=self.genVars_From_Term(Term)
        f=self.transTerm_to_ANF(Term)
        corr_list=[0 for i in range(0,2**len(Cut_points))]
        if set(Vars) & set(Cut_points)!=set():
            for cut_len in range(0,2**len(Cut_points)):
                #cut_values=list(bin(cut_len)[2:].zfill(len(Cut_points)))
                cut_values=list(self.wordToBinaryString(cut_len, len(Cut_points)))
                f_tem=f
                for k in range(0,len(Cut_points)):
                    f_tem=f_tem.subs([(symbols(Cut_points[k]),int(cut_values[k]))])
                corr_list[cut_len]=self.correlation_of_smallF(f_tem, list(set(Vars)-set(Cut_points))) 
            return corr_list
        else:
            corr=self.correlation_of_smallF(f, Vars)
            return [corr]
        
    
    def justify_Term(self,Term):
        """
        To ensure graph have no single node
        """
        G=self.genGraph(Term)
        G_subs=list(nx.connected_components(G))
        if min([len(list(c)) for c in G_subs])==1:
            return 0
        else:
            return 1
        
    def correlation_of_Term(self,Term,file):
        t1=time.time()
        global flag
        myfile=open(file,'a')
        Term=self.reduce_Term(Term)
        
        # If Graph has single node, correlation is 0
        if self.justify_Term(Term)==0:
            print("correlation is: 0",file=myfile)
            print("----------",file=myfile)
            t2=time.time()
            print("Time used:"+str(t2-t1)+"s",file=myfile)
            print(file=myfile)
            return 0
        
        else:
            Vars=self.genVars_From_Term(Term)
            if len(Vars)<=self.Threshold:
                print("correlation is:",self.correlation_of_smallTerm(Term),file=myfile)
                return self.correlation_of_smallTerm(Term)
            else:
                Cut_points,divided_Term=self.divide_Term(Term)
                print('Cuts:',len(Cut_points),Cut_points)
                for v in divided_Term:
                    print(len(v),v)
                print("Cut_points are:",Cut_points,file=myfile)
                print("----------",file=myfile)
                cut_Corr=[1 for i in range(0,2**len(Cut_points))]
                corr=1
                print("divided_Terms are:",file=myfile)
                for i in range(0,len(divided_Term)):
                    Vars=self.genVars_From_Term(divided_Term[i])
                    corr_list=self.correlation_of_Term_with_cutpoints(divided_Term[i], Cut_points,myfile)
                    if len(corr_list)==1:
                        corr=corr*corr_list[0]
                    elif len(corr_list)>1:
                        cut_Corr = list(map(lambda x: x[0]*x[1], zip(cut_Corr, corr_list)))  # Piling-up lemma
                    if ((abs(max(cut_Corr))+abs(min(cut_Corr)))==0)|(corr==0):
                        print("----------",file=myfile)
                        print("correlation is: 0",file=myfile)
                        t2=time.time()
                        print("Time used:"+str(t2-t1)+"s",file=myfile)
                        print(file=myfile)
                        return 0
                print("----------",file=myfile)
                cut_Corr=list(map(lambda x: x*corr,cut_Corr))
                corr=np.mean(cut_Corr)* (-1)**flag
                print("flag is: ",flag, file=myfile)
                print("sub-correlation are:",file=myfile)
                print("CutPoints' value",file=myfile)
                for cut_len in range(0,2**len(Cut_points)):
                    cut_values=bin(cut_len)[2:].zfill(len(Cut_points))
                    if cut_Corr[cut_len]!=0.0:
                        print(cut_values,"\t\t\t"+str((-1)**(int(cut_Corr[cut_len]<0))*2)+"^",math.log2(abs(cut_Corr[cut_len])),\
                              file=myfile)
                    elif cut_Corr[cut_len]==0:
                        print(cut_values,'\t\t\t','0',file=myfile)
                print("----------",file=myfile)
                if corr>0:
                    print("correlation is: 2^",math.log2(abs(corr)),\
                      '(',corr,')',file=myfile)
                elif corr<0:
                    print("correlation is: - 2^",math.log2(abs(corr)),\
                      '(',corr,')',file=myfile)
                else:
                    print("correlation is: 0",file=myfile)
                print("----------",file=myfile)
                t2=time.time()
                print("Time used:"+str(t2-t1)+"s",file=myfile)
                print(file=myfile)
                return corr

        