
import numpy as np
import copy

"""Define global variables"""
computed_table = dict()
inique_table = dict()
global_var_order = []
qubit_num = 0

class Node:
    """The node of QMDD"""
    def __init__(self,key):
        self.key = key
        self.idx=0
        self.successor=[None]*4
        
class QMDD:
    """The node of QMDD"""
    def __init__(self,node):
        self.weight=1
        if isinstance(node,Node):
            self.node=node
        else:
            self.node=Node(node)
    
    def node_number(self):
        node_set=set()
        node_set=get_node_set(self,node_set)
        return len(node_set)
    
    def self_copy(self):
        temp=QMDD(self.node)
        temp.weight=self.weight
        return temp
    
    def __eq__(self,other):
        if self.node==other.node and get_int_key(self.weight)==get_int_key(other.weight):
            return True
        else:
            return False
    
def get_node_set(qmdd,node_set=set()):
    """Only been used when counting the node number of a QMDD"""
    node_set.add(qmdd.node)
    for k in range(4):
        if qmdd.node.successor[k]:
            node_set=get_node_set(qmdd.node.successor[k],node_set)
    return node_set
    

def get_int_key(weight):
    """To transform a complex number to a tuple with int values"""
    epi=0.00001     
    return (int(round(weight.real/epi)) ,int(round(weight.imag/epi)))

def find_computed_table(qmdd1,qmdd2,op):
    """To return the results that already exist"""
    global computed_table
    the_key = (get_int_key(qmdd1.weight),qmdd1.node,get_int_key(qmdd2.weight),qmdd2.node,op)
    if computed_table.__contains__(the_key):
        res = computed_table[the_key]
        qmdd = QMDD(res[1])
        qmdd.weight = res[0]
        return qmdd
    else:
        return None

def insert_2_computed_table(qmdd1,qmdd2,op,res):
    """To insert an item to the computed table"""
    global computed_table
    the_key = (get_int_key(qmdd1.weight),qmdd1.node,get_int_key(qmdd2.weight),qmdd2.node,op)
    computed_table[the_key] = (res.weight,res.node)

def Find_Or_Add_Unique_table(x,succ):
    """To return a node if it is already exist, create a new node otherwise"""    
    global unique_table
    
    if not isinstance(x,str):
        if unique_table.__contains__(x):
            return unique_table[x]
        else:
            res=Node(x)
            unique_table[x]=res
        return res
    
    weig_in=[get_int_key(c.weight) for c in succ]

    temp_key=(x,weig_in[0],weig_in[1],weig_in[2],weig_in[3],succ[0].node,succ[1].node,succ[2].node,succ[3].node)
    if unique_table.__contains__(temp_key):
        return unique_table[temp_key]
    else:
        res=Node(x)
        res.successor = succ
        unique_table[temp_key]=res
    return res 

def normalize(x,succ):
    """The normalize and reduce procedure"""
    weig_in=[get_int_key(c.weight) for c in succ]
    
    if weig_in[0]==weig_in[1]==weig_in[2]==weig_in[3]==(0,0):
        temp=QMDD(Find_Or_Add_Unique_table(1,[]))
        temp.weight=0
        return temp
    for k in range(4):
        if weig_in[k]==(0,0):
            succ[k]=QMDD(Find_Or_Add_Unique_table(1,[]))
            succ[k].weight=0
      
    epi=0.00001        
    weig_abs=[np.around(abs(c.weight)/epi) for c in succ]
    if max(weig_abs)!=0:
        w = succ[weig_abs.index(max(weig_abs))].weight
    else:
        w=1
       
    for k in range(4):
        succ[k].weight/=w
    node = Find_Or_Add_Unique_table(x,succ)
    res = QMDD(node)
    res.weight = w
    return res


def Ini_QMDD(n):
    """Initialize the package with qubit number n;
    """

    global unique_table    
    unique_table=dict()
    global computed_table
    computed_table = dict()
    global global_var_order
    order=[]
    for k in range(n-1,-1,-1):
        order.append('x'+str(k))
        
    global_var_order = order
    global qubit_num
    qubit_num = n
    
    I=np.array([[1,0],
           [0,1]])
    
    qmdd=get_single_QMDD(I,0)
    return qmdd    
    

def get_single_QMDD(U,qubit):
    """To get the QMDD of a single qubit gate"""
    I=np.array([[1,0],
            [0,1]])
    
    global qubit_num
    n=qubit_num
    
    mid_res=matrix_QMDD(U,qubit)
    if not isinstance(qubit,int): 
        q=qubit[0] 
    else:
        q=qubit
    for k in range(q+1,n,1):
        temp=matrix_QMDD(I,k)
        mid_res=QMDD_Tensor(temp,mid_res)
    for k in range(q-1,-1,-1):
        temp=matrix_QMDD(I,k)
        mid_res=QMDD_Tensor(mid_res,temp)

    return mid_res

def get_cnot_QMDD(q):
    """To get the QMDD of a CNOT gate but not expand to a n-qubit form， q=[ctrl, targ]"""
    
    term=Find_Or_Add_Unique_table(1,[])
        
    I=np.array([[1,0],
               [0,1]])
    X=np.array([[0,1],
               [1,0]])
    if q[0]>q[1]:
        x='x'+str(q[0])
        succ=[None]*4
        succ[1]=QMDD(term)
        succ[1].weight=0
        succ[2]=QMDD(term)
        succ[2].weight=0
        mid_res=matrix_QMDD(I,q[1])
        for k in range(q[1]+1,q[0],1):
            temp=matrix_QMDD(I,k)
            mid_res=QMDD_Tensor(temp,mid_res)
        succ[0]=mid_res
 
        mid_res2=matrix_QMDD(X,q[1])
        for k in range(q[1]+1,q[0],1):
            temp=matrix_QMDD(I,k)
            mid_res2=QMDD_Tensor(temp,mid_res2)
        succ[3]=mid_res2

        
    if q[0]<q[1]:
        x='x'+str(q[1])
        succ=[None,None,None,None]
        temp_succ=[QMDD(term),QMDD(term),QMDD(term),QMDD(term)]
        temp_succ[1].weight=0
        temp_succ[2].weight=0
        temp_succ[3].weight=0
        node = Find_Or_Add_Unique_table('x'+str(q[0]),temp_succ)
        temp1=QMDD(node)
        for k in range(q[0]+1,q[1],1):       
            temp=matrix_QMDD(I,k)
            temp1=QMDD_Tensor(temp,temp1)
        succ[0]=temp1
        succ[3]=temp1.self_copy()
        temp_succ=[QMDD(term),QMDD(term),QMDD(term),QMDD(term)]
        temp_succ[0].weight=0
        temp_succ[1].weight=0
        temp_succ[2].weight=0
        node = Find_Or_Add_Unique_table('x'+str(q[0]),temp_succ)
        temp2=QMDD(node)
        for k in range(q[0]+1,q[1],1):       
            temp=matrix_QMDD(I,k)
            temp2=QMDD_Tensor(temp,temp2)
        succ[1]=temp2
        succ[2]=temp2.self_copy()
        
    return normalize(x,succ)
    
    
    
    
def get_cnot_QMDD_n(q):
    """To get the QMDD of a CNOT gate and expand it a n-qubit form， q=[ctrl, targ]"""
    global qubit_num
    temp=get_cnot_QMDD(q)
    
    n=qubit_num
    
    I=np.array([[1,0],
             [0,1]])
        
    for k in range(max(q)+1,n,1):
        temp_I=matrix_QMDD(I,k)
        temp=QMDD_Tensor(temp_I,temp)
    for k in range(min(q)-1,-1,-1):
        temp_I=matrix_QMDD(I,k)
        temp=QMDD_Tensor(temp,temp_I)
    return temp
    
    
    
def matrix_QMDD(U,q):
    """To get the QMDD of a matrix U on the qubit list q"""
        
    (m,n)=U.shape
    
    if m==n==1:
        node=Find_Or_Add_Unique_table(1,[])
        temp=QMDD(node)
        temp.weight=U[0][0]
        return temp
    
    succ=[None]*4

    if isinstance(q,int):
        if m==n==2:
            x='x'+str(q)
            for k in range(4):
                mid_res=matrix_QMDD(U[k//2:k//2+1,k%2:k%2+1],q)
                succ[k]=mid_res
        else:
            raise(Exception('Not the correct matrix size'))
    else:
        x='x'+str(q[0])
        for k in range(4):
            mid_res=matrix_QMDD(U[(k//2)*m//2:(k//2+1)*m//2,(k%2)*n//2:(k%2+1)*n//2],q[1:])
            succ[k]=mid_res

    return normalize(x,succ)      


def QMDD_Tensor(qmdd1,qmdd2):
    """The Tensor product of qmdd1 and qmdd2"""
    w1=qmdd1.weight
    v1=qmdd1.node.key
    w2=qmdd2.weight
    v2=qmdd2.node.key  
    if find_computed_table(qmdd1,qmdd2,'tensor'):
        return find_computed_table(qmdd1,qmdd2,'tensor')
    
    if not isinstance(v1,str):
        if w1==0:
            return qmdd1.self_copy()
        if w1==1:
            return qmdd2.self_copy()
        temp=QMDD(qmdd2.node)
        temp.weight=w1*w2
        return temp
    succ=[None]*4

    for k in range(4):
        p=qmdd1.node.successor[k]
        q=qmdd2
        mid_res=QMDD_Tensor(p,q)
        succ[k]=mid_res
    temp=normalize(v1,succ)
    temp.weight*=w1
    insert_2_computed_table(qmdd1,qmdd2,'tensor',temp)
    return temp


def QMDD_ADD(qmdd1,qmdd2):
    """The addtition of two QMDDs;
    return qmdd1+qmdd2
    """
    
    if find_computed_table(qmdd1,qmdd2,'+'):
        return find_computed_table(qmdd1,qmdd2,'+')
    
    global global_var_order
    
    w1=qmdd1.weight
    v1=qmdd1.node.key
    w2=qmdd2.weight
    v2=qmdd2.node.key
    
    
    if not isinstance(v1,str):
        if w1==0:
            return qmdd2.self_copy()
        if not isinstance(v2,str):
            temp=QMDD(qmdd1.node)
            temp.weight=w1+w2
            insert_2_computed_table(qmdd1,qmdd2,'+',temp)
            return temp
        idx_v1=float('inf')
    else:
        idx_v1=global_var_order.index(v1)
        
    if not isinstance(v2,str):
        if w2==0:
            return qmdd1.self_copy()
        idx_v2=float('inf')
    else:
        idx_v2=global_var_order.index(v2)
        
    succ=[None]*4

    if v1==v2:
        x=v1
        for k in range(4):
            p=qmdd1.node.successor[k].self_copy()
            p.weight*=w1
            q=qmdd2.node.successor[k].self_copy()
            q.weight*=w2
            mid_res=QMDD_ADD(p,q)
            succ[k]=mid_res
            
    if idx_v1<idx_v2:
        x=v1
        for k in range(4):
            p=qmdd1.node.successor[k].self_copy()
            p.weight*=w1
            q=qmdd2
            mid_res=QMDD_ADD(p,q)
            succ[k]=mid_res
            
    if idx_v1>idx_v2:
        x=v2
        for k in range(4):
            p=qmdd1
            q=qmdd2.node.successor[k].self_copy()
            q.weight*=w2
            mid_res=QMDD_ADD(p,q)
            succ[k]=mid_res
            
    temp=normalize(x,succ)
    insert_2_computed_table(qmdd1,qmdd2,'+',temp)
    return temp   


def QMDD_Mul(qmdd1,qmdd2):
    """The multiplication of two QMDDs;
    return qmdd1*qmdd2
    
    """
    if find_computed_table(qmdd1,qmdd2,'*'):
        return find_computed_table(qmdd1,qmdd2,'*')    
    
    global global_var_order
    w1=qmdd1.weight
    v1=qmdd1.node.key
    w2=qmdd2.weight
    v2=qmdd2.node.key
    
    if not isinstance(v1,str):
        if w1==0:
            return qmdd1.self_copy()
        temp=QMDD(qmdd2.node)
        temp.weight=w1*w2
        insert_2_computed_table(qmdd1,qmdd2,'*',temp)
        return temp
        
    if not isinstance(v2,str):
        if w2==0:
            return qmdd2.self_copy()
        temp=QMDD(qmdd1.node)
        temp.weight=w1*w2
        insert_2_computed_table(qmdd1,qmdd2,'*',temp)
        return temp

    succ=[None]*4
    if v1==v2:
        x=v1
        for k1 in range(2):
            for k2 in range(2):
                mid_res=QMDD(Find_Or_Add_Unique_table(1,[]))
                mid_res.weight=0
                for k3 in range(2):
                    p=qmdd1.node.successor[2*k1+k3]
                    q=qmdd2.node.successor[k2+2*k3]
                    mid_res=QMDD_ADD(QMDD_Mul(p,q),mid_res)
                succ[2*k1+k2]=mid_res
        temp=normalize(x,succ)
        temp.weight*=w1*w2
        insert_2_computed_table(qmdd1,qmdd2,'*',temp)
        return temp                
                
    if global_var_order.index(v1)<global_var_order.index(v2):
        x=v1
        for k1 in range(2):
            for k2 in range(2):
                mid_res=QMDD(Find_Or_Add_Unique_table(1,[]))
                mid_res.weight=0
                for k3 in range(2):
                    p=qmdd1.node.successor[2*k1+k3]
                    q=qmdd2
                    mid_res=QMDD_ADD(QMDD_Mul(p,q),mid_res)
                succ[2*k1+k2]=mid_res
        temp=normalize(x,succ)
        temp.weight*=w1
        insert_2_computed_table(qmdd1,qmdd2,'*',temp)
        return temp                

                
    if global_var_order.index(v1)>global_var_order.index(v2):
        x=v2
        for k1 in range(2):
            for k2 in range(2):
                mid_res=QMDD(Find_Or_Add_Unique_table(1,[]))
                mid_res.weight=0
                for k3 in range(2):
                    p=qmdd1
                    q=qmdd2.node.successor[k2+2*k3]
                    mid_res=QMDD_ADD(QMDD_Mul(p,q),mid_res)
                succ[2*k1+k2]=mid_res
        temp=normalize(x,succ)
        temp.weight*=w2
        insert_2_computed_table(qmdd1,qmdd2,'*',temp)
        return temp                