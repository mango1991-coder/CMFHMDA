import xlrd
import numpy
import numpy as np                                
import math                                              
import heapq                                       
from datetime import datetime
start0 = datetime.now()
import pandas as pd
import numpy.linalg as LA

file3 = 'disease-microbe associations.xlsx'

a = open(r'disease-microbe associations.xlsx')
b = open(r'Disease symptom similarity.xlsx')
A=numpy.zeros((39,292))
zhengzhuangleisixing = numpy.zeros((39, 39))
xlsx1 = xlrd.open_workbook(r'disease-microbe associations.xlsx')
sheet1 = xlsx1.sheets()[0]
for i in range(450):
    s1=sheet1.row_values(i)
    m=int(s1[0])
    n=int(s1[1])
    A[m-1,n-1]=1  
xlsx2 = xlrd.open_workbook(r'Disease symptom similarity.xlsx')
sheet2 = xlsx2.sheets()[0]
for i in range(61):
    s2 = sheet2.row_values(i)
    z1 = int(s2[0])
    z2 = int(s2[1])
    z3 = float(s2[2])
    zhengzhuangleisixing[z1 - 1, z2 - 1] = zhengzhuangleisixing[z2 - 1, z1 - 1] = z3
for i in range(39):
    zhengzhuangleisixing[i,i]=1

### S_M
table3 = xlrd.open_workbook(filename=file3).sheets()[0]    
row3 = table3.nrows  
col3 = table3.ncols  
S_M = np.zeros((row3, col3))                 
for z in range(col3):
    cols = np.matrix(table3.col_values(z))    
    S_M[:,z] = cols                          
    
Y0=np.zeros((39,292))
for i in range(450):
    Y0[int(S_M[i,0])-1,int(S_M[i,1])-1]=1 

M = len(Y0)
N = len(Y0[0]) 
####  
k = 5                                           
K = 39                                             
P = 0.7                                           
lamda1 = 2
lamdam = 0.1
lamdad = 0.05
Ik = np.identity(K)                                
##########   WKNKN                 
def WKNKN(Y0,S_S,M_M):                                             
    distances=np.zeros((M,M))
    distancem=np.zeros((N,N))
    global Ys,Ym
    Ys=np.zeros((M,N))   
    Ym=np.zeros((M,N))         
    for i in range(M):
        for j in range(M):                                                               
            if j!=i: 
                distances[i][j]=np.linalg.norm(S_S[i,:]-S_S[j,:])                  
        distance_s=distances[i,:].tolist()                                         
        p=list(map(distance_s.index,heapq.nsmallest(k+1,distance_s)))            
        SS_num=0
        Ws=0
        for a in range(1,k+1):
            Ws +=math.pow(P,a-1)*S_S[i][p[a]]*Y0[p[a],:]                            
            SS_num+=S_S[i][p[a]]                                              
        Ys[i,:]=Ws/SS_num                                                             
#########################  
    for j in range(N):
        for i in range(N):                          
            if i!=j:              
                distancem[i][j]=np.linalg.norm(M_M[:,i]-M_M[:,j])                 
        distance_m=distancem[:,j].tolist()                                         
        q=list(map(distance_m.index,heapq.nsmallest(k+1,distance_m)))           
        SM_num=0
        Wm=0                                                                       
        for b in range(1,k+1):                            
            Wm +=math.pow(P,b-1)*M_M[q[b]][j]*Y0[:,q[b]]                            
            SM_num+=M_M[q[b]][j]
        if j !=291:
            Ym[:,j]=Wm/SM_num                                                       
    global Ysm        
    Ysm=(Ys+Ym)/2
    Y_WKNKN=np.maximum(Y0,Ysm)   
    return Y_WKNKN        
##############################   WKNKN  

for i_new in range(450):
    print ('rank',i_new+1,'times')
    start = datetime.now()
    score_result_0=[]
    
    Y=Y0.copy()
    Y[int(S_M[i_new,0]-1),int(S_M[i_new,1]-1)]=0
    C = np.asmatrix(Y)
    gamd = 39 / (LA.norm(C, 'fro') ** 2);
    kd = np.mat(np.zeros((39, 39)))
    km = np.mat(np.zeros((292, 292)))
    D = C * C.T;
    for i in range(39):
        for j in range(i, 39):
            kd[j, i] = np.exp(-gamd * (D[i, i] + D[j, j] - 2 * D[i, j]))
    kd = kd + kd.T - np.diag(np.diag(kd))
    KD = np.asarray(kd)
    kd = []  
    SD = (KD + zhengzhuangleisixing) / 2
    S_S = np.asarray(SD)
    df = pd.DataFrame(SD)
    excel_filename = 'disease similarity matrix.xlsx'
    df.to_excel(excel_filename, index=False, header=False)
    gamam = 292 / (LA.norm(C, 'fro') ** 2);
    E = C.T * C;
    for i in range(292):
        for j in range(i, 292):
            km[i, j] = np.exp(-gamam * (E[i, i] + E[j, j] - 2 * E[i, j]))
    #    kd=kd+kd.T-np.diag(np.diag(kd))
    km = km + km.T - np.diag(np.diag(km))
    M_M = np.asarray(km)
    df = pd.DataFrame(km)
    excel_filename1 = 'microbe similarity matrix.xlsx'
    df.to_excel(excel_filename1, index=False, header=False)
    km = []
    file1 = 'microbe similarity matrix.xlsx'
    file2 = 'disease similarity matrix.xlsx'
    ### M_M
    table1 = xlrd.open_workbook(filename=file1).sheets()[0]
    row1 = table1.nrows
    col1 = table1.ncols
    M_M = np.zeros((row1, col1))
    for x in range(row1):
        rows = np.matrix(table1.row_values(x))
        M_M[x, :] = rows
        ### S_S
    table2 = xlrd.open_workbook(filename=file2).sheets()[0]
    row2 = table2.nrows
    col2 = table2.ncols
    S_S = np.zeros((row2, col2))
    for y in range(col2):
        cols = np.matrix(table2.col_values(y))
        S_S[:, y] = cols

    #Y_WKNKN=WKNKN(Y,S_S,M_M)
    Y_WKNKN=Y

                 
    ##########    SVD
    U,Sk,V = np.linalg.svd(np.matrix(Y_WKNKN), full_matrices=False)
    ###U, Sk, V = np.linalg.svd(np.matrix(Y), full_matrices=False)
    S = Sk[:K]                                                     
    
    UK = U[:,0:K]                                                                                                                  
    VK = V[0:K,:]                                                  
    
    sigma = np.diag(S)                                             
    Sigm = np.matrix(sigma)                                        
    Sigma = np.sqrt(Sigm)                                       
            
    VK = VK.T                                                      
    A = UK*Sigma
    B = VK*Sigma                                        
    ########################## 
    for T0 in range(1000):                                                  
        A = (Y_WKNKN*B + lamdad*S_S*A)*np.linalg.inv(B.T*B+lamda1*Ik+lamdad*A.T*A)
        B = (Y_WKNKN.T*A + lamdam*M_M*B)*np.linalg.inv(A.T*A+lamda1*Ik+lamdam*B.T*B)
    ########################## 
    
    Y_hat = np.dot(A, B.T)

    df_Y = pd.DataFrame(Y_hat)

    df_Y.to_excel('1(lamdad = 0.05).xlsx', index=False, engine='openpyxl')


    for j in range(N):
        if Y0[int(S_M[i_new,0])-1,j]==0:
            score_result_0.append(Y_hat[int(S_M[i_new,0])-1,j])
    score_result = Y_hat[int(S_M[i_new,0])-1,int(S_M[i_new,1])-1]  
    score_result_0.append(score_result)          
    score2 = sorted(score_result_0,reverse=True)                                          
    rank_result = np.average(np.where(np.array(score2)==score_result))+1      
               
###                      
    print('result=%d,   local2=%d'%(rank_result))
    print('\n')    

    
    result=open('result(lamdad = 0.05).txt','a+')
    result.writelines([str(rank_result),'\t'])
    result.close()
           
    end = datetime.now()
    print((end-start))

end0 = datetime.now()
print((end0-start0))








