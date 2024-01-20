'''
Created on 20190830

@author: Administrator
'''
import sys
sys.path.append('..')

from ToolsGrain import *

A = [1,2,4,10,31,43,56]

def genTerms_Tz_Tb(Tz,Tb):
    
    Terms = []

    # g'(b^(t+j))
    for j in A:
        for i in Tz:
            Term = [['b_' + str(j + i)]]
            Terms = Terms + Term
    
    # b_{t+128+j}, g(b^(t+j))
    for j in Tb:
        Term = [['b_' + str(j + 80)], \
          ['b_' + str(j + 62)], ['b_' + str(j + 60)], ['b_' + str(j + 52)], ['b_' + str(j + 45)], ['b_' + str(j + 37)], \
          ['b_' + str(j + 33)],['b_' + str(j + 28)],['b_' + str(j + 21)],['b_' + str(j + 14)],['b_' + str(j + 9)],\
          ['b_' + str(j + 0)],['b_' + str(j + 63), 'b_' + str(j + 60)], ['b_' + str(j + 37), 'b_' + str(j + 33)],\
          ['b_' + str(j + 15), 'b_' + str(j + 9)],['b_' + str(j + 60), 'b_' + str(j + 52),'b_' + str(j + 45)],\
          ['b_' + str(j + 33),'b_' + str(j + 28),'b_' + str(j + 21)],\
          ['b_' + str(j + 63),'b_' + str(j + 45),'b_' + str(j +28),'b_' + str(j + 9)],\
          ['b_' + str(j + 60),'b_' + str(j + 52),'b_' + str(j + 37),'b_' + str(j + 33)],\
          ['b_' + str(j + 63),'b_' + str(j + 60),'b_' + str(j + 21),'b_' + str(j + 15)],\
          ['b_' + str(j + 63),'b_' + str(j + 60),'b_' + str(j + 52),'b_' + str(j + 45),'b_' + str(j + 37)],\
          ['b_' + str(j + 33),'b_' + str(j + 28),'b_' + str(j + 21),'b_' + str(j + 15),'b_' + str(j + 9)],\
          ['b_' + str(j + 52),'b_' + str(j + 45),'b_' + str(j + 37),'b_' + str(j + 33),'b_' + str(j + 28),'b_' + str(j + 21)]]
        Terms = Terms + Term
    return Terms

Tz = [0,9,14,21,28,33,37,45,52,60,62,80]  # corr2
Tb = [1,2,4,10,31,43,56] # corr2

# Tz=[0,14,21,60,62,80] # corr1
# Tb=[1,2,4,10,31,43,56] # corr1

# Tz=[0,14,21,28,37,45,52,60,62,80] # corr3
# Tb=[1,2,4,10,31,43,56] # corr3

Terms = genTerms_Tz_Tb(Tz, Tb)


Terms_Tz = []
Gamma_name = []
for i in Tz:
    Terms_Tz = Terms_Tz + [['b_' + str(i + 63)] ]
    Gamma_name = Gamma_name + ['G_' + str(i) + '[' + str(0) + ']'] 
for v in Terms_Tz:
    print('Terms:',v)

time1=time.time()
T = Correlation_of_Term_F(10, 1)
# Terms=Terms+Terms_Tz
Terms = T.reduce_Term(Terms)


for v in Terms:print(v)
G = T.genGraph(Terms)
T.showGraph(G)
fileProcess = "CG" + ".txt"
# corr = T.correlation_of_Term(Terms, fileProcess)

  
G_cutsets,divided_Term=T.divide_Term(Terms)
print('-----------Cut Sets-----------')
print(len(G_cutsets),G_cutsets)
print('-----------End Cut Sets-----------')
print('-----------Divided Terms Corr-----------')
for v in divided_Term:
    print(len(v),v)

Corr_list=[[] for i in range(0,len(divided_Term))]
divided_vars=[]
for i in range(0,len(divided_Term)):
    print('divided_Term',i)
    Corr_list[i]=T.correlation_of_Term_with_cutpointsNonFile(divided_Term[i],G_cutsets)
    divided_vars=divided_vars+[T.genVars_From_Term(divided_Term[i])]
print('-----------End Divided Terms Corr-----------')
     

#Term_H=[['b_77']] # corr1
Term_H=[['b_77'],['b_84'],['b_91'],['b_108']]  # corr3
H_flag=[None for i in range(0,len(Term_H))]
for i in range(0,len(Term_H)):
    for j in range(0,len(divided_vars)):
        if Term_H[i][0] in divided_vars[j]:
            H_flag[i]=j

G_corr_Value=[]
for i in range(0,2**len(Term_H)):
    Corr_list_temp=Corr_list+[]
    Corr_end=[1 for i in range(0,2**len(G_cutsets))]
    term_temp=[[] for j in range(0,len(divided_Term))]
    H=list(Basic_Tools.wordToBinaryString(i, len(Term_H)))
    print(H)
    for j in range(0,len(Term_H)):
        if H[j] == '1':
            term_temp[H_flag[j]]=term_temp[H_flag[j]]+[Term_H[j]]
    for j in range(0,len(divided_Term)):
        if term_temp[j]!=[]:
            Corr_list_temp[j]=T.correlation_of_Term_with_cutpointsNonFile((divided_Term[j]+term_temp[j]),G_cutsets)
        Corr_end = list(map(lambda x: x[0]*x[1], zip(Corr_end, Corr_list_temp[j])))
    corr=np.mean(Corr_end)
    G_corr_Value=G_corr_Value+[corr]
    if corr<0:
        print(" - 2^",math.log2(abs(corr)))
    if corr>0:
        print(" 2^",math.log2(abs(corr))) 
    if corr==0:
        print(0)    
time2=time.time()
print('Time used: '+str(time2-time1)+'s')
        
print(G_corr_Value)  



'''
# corr2
Term_H=[['b_77'],['b_84'],['b_91'],['b_108']]  # corr2
H_flag=[None for i in range(0,4)]
for i in range(0,4):
    for j in range(0,len(divided_vars)):
        if Term_H[i][0] in divided_vars[j]:
            H_flag[i]=j

G_corr_Value=[]
for i in range(0,2**4):
    Corr_list_temp=Corr_list+[]
    Corr_end=[1 for i in range(0,2**len(G_cutsets))]
    term_temp=[[] for j in range(0,len(divided_Term))]
    H=list(Basic_Tools.wordToBinaryString(i, 4))
    print(H)
    for j in range(0,4):
        if H[j] == '1':
            term_temp[H_flag[j]]=term_temp[H_flag[j]]+[Term_H[j]]
    for j in range(0,len(divided_Term)):
        if term_temp[j]!=[]:
            Corr_list_temp[j]=T.correlation_of_Term_with_cutpointsNonFile((divided_Term[j]+term_temp[j]),G_cutsets)
        Corr_end = list(map(lambda x: x[0]*x[1], zip(Corr_end, Corr_list_temp[j])))
    corr=np.mean(Corr_end)
    G_corr_Value=G_corr_Value+[corr]
    if corr<0:
        print(" - 2^",math.log2(abs(corr)))
    if corr>0:
        print(" 2^",math.log2(abs(corr))) 
    if corr==0:
        print(0)    
time2=time.time()
print('Time used: '+str(time2-time1)+'s')
        
print(G_corr_Value)  

G_corr_Index=[]
for i in range(0,2**4):
    H=list(Basic_Tools.wordToBinaryString(i, 4))
    G_corr_Index=G_corr_Index+[[0,0,int(H[0]),int(H[1]),int(H[2]),0,0,int(H[3]),0,0,0,0]]
for v in G_corr_Index:
    print()
'''


