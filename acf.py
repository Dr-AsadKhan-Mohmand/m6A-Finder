from __future__ import division
import numpy as np
from statsmodels.tsa.stattools import acf
from Bio import SeqIO
#from sklearn.ensemble import GradientBoostingClassifier
#from sklearn.model_selection import cross_val_score
#from sklearn.model_selection import KFold
#from sklearn.svm import SVC
nl=15
nparr = np.array([])
for rec in SeqIO.parse("m6.fasta", "fasta"):
    str1="".join(rec.seq)
    A=[0.2837,0.3183,0.1873,0.3299,0.2652,0.3696,0.1292,0.2618]
    C=[0.2162,0.2098,0.2602,0.2031,0.2180,0.1901,0.2902,0.2319]
    G=[0.2972,0.2643,0.2882,0.2668,0.2966,0.2661,0.2572,0.2810]
    U=[0.2027,0.2073,0.2642,0.2000,0.2200,0.1741,0.3231,0.2251]
    S=[]
    P1=[]
    P2=[]
    P3=[]
    P4=[]
    P5=[]
    P6=[]
    P7=[]
    P8=[]
    a=0
    c=0
    g=0
    u=0
    for i in range (0,len(str1)):
        if(str1[i]=="A"):
            P1.append(A[0])
            P2.append(A[1])
            P3.append(A[2])
            P4.append(A[3])
            P5.append(A[4])
            P6.append(A[5])
            P7.append(A[6])
            P8.append(A[7])
        if(str1[i]=="C"):
            P1.append(C[0])
            P2.append(C[1])
            P3.append(C[2])
            P4.append(C[3])
            P5.append(C[4])
            P6.append(C[5])
            P7.append(C[6])
            P8.append(C[7])
        if(str1[i]=="G"):
            P1.append(G[0])
            P2.append(G[1])
            P3.append(G[2])
            P4.append(G[3])
            P5.append(G[4])
            P6.append(G[5])
            P7.append(G[6])
            P8.append(G[7])
        if(str1[i]=="U"):
            P1.append(U[0])
            P2.append(U[1])
            P3.append(U[2])
            P4.append(U[3])
            P5.append(U[4])
            P6.append(U[5])
            P7.append(U[6])
            P8.append(U[7])
    a1=acf(P1,nlags=nl)
    a2=acf(P2,nlags=nl)
    a3=acf(P3,nlags=nl)
    a4=acf(P4,nlags=nl)
    a5=acf(P5,nlags=nl)
    a6=acf(P6,nlags=nl)
    a7=acf(P7,nlags=nl)
    a8=acf(P8,nlags=nl)
    Sr=(a1+a2+a3+a4+a5+a6+a7+a8)/8
    nparr = np.append(nparr, np.array(Sr),axis=0)
    #print(len(nparr))
    #if(len(nparr)%204==0):
        #print("true")

resh=nparr.reshape(2614,nl+1)
matrix=np.delete(resh,0,1)
np.savetxt("C:acf.csv",matrix, delimiter=",")

        