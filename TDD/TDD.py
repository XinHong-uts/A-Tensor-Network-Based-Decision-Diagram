import numpy as np
import copy

"""Define global variables"""
computed_table = dict()
unique_table = dict()
global_index_order = []


class Index:
    """The index, here idx is used when there is a hyperedge"""
    def __init__(self,key,idx=0):
        self.key = key
        self.idx = idx
        
    def __eq__(self,other):
        if self.key == other.key and self.idx == other.idx:
            return True
        else:
            return False
    def __lt__(self,other):
        self_order = global_index_order.index(self.key)
        other_order = global_index_order.index(other.key)
        return self_order < other_order
    
class Node:
    """To define the node of TDD"""
    def __init__(self,key):
        self.key = key
        self.idx = 0
        self.out_weight=[1]*2
        self.successor=[None]*2

class TDD:
    def __init__(self,node):
        """TDD"""
        self.weight=1
        
        self.index_set=[]
        
        if isinstance(node,Node):
            self.node=node
        else:
            self.node=Node(node)
            
    def node_number(self):
        node_set=set()
        node_set=get_node_set(self.node,node_set)
        return len(node_set)
    
    def self_copy(self):
        temp = TDD(self.node)
        temp.weight = self.weight
        temp.index_set = self.index_set
        return temp
    
    def __eq__(self,other):
        if self.node==other.node and get_int_key(self.weight)==get_int_key(other.weight):
            return True
        else:
            return False



def Ini_TDD(var_order):
    """To initialize the unique_table,computed_table and set up a global index order"""
    global computed_table
    global unique_table
    global global_index_order
    
    unique_table = dict()
    computed_table = dict()
    global_index_order = var_order
    
    return True

def get_int_key(weight):
    """To transform a complex number to a tuple with int values"""
    epi=0.00001     
    return (int(round(weight.real/epi)) ,int(round(weight.imag/epi)))

def get_node_set(node,node_set=set()):
    """Only been used when counting the node number of a TDD"""
    node_set.add(node)
    for k in range(2):
        if node.successor[k]:
            node_set = get_node_set(node.successor[k],node_set)
    return node_set

def normalize(x,low,high):
    """The normalize and reduce procedure"""
    if high==low:
        return high
    weig1=low.weight
    weig2=high.weight
    if get_int_key(weig1)==(0,0) and get_int_key(weig2)==(0,0):
        node=Find_Or_Add_Unique_table(1,0,0,None,None)
        res=TDD(node)
        res.weight=0
        return res
    elif get_int_key(weig1)==(0,0):
        node=Find_Or_Add_Unique_table(1,0,0,None,None)
        low=TDD(node)
        weig1=0
    elif get_int_key(weig2)==(0,0):
        node=Find_Or_Add_Unique_table(1,0,0,None,None)
        high=TDD(node)
        weig2=0
        
    epi=0.00001
    if np.around(abs(weig1)/epi)>=np.around(abs(weig2)/epi):
        weig=weig1
    else:
        weig=weig2

    weig1=weig1/weig
    weig2=weig2/weig
    node=Find_Or_Add_Unique_table(x,weig1,weig2,low.node,high.node)
    res=TDD(node)
    res.weight=weig
    return res

def Find_Or_Add_Unique_table(x,weight1,weight2,node1,node2):
    """To return a node if it already exist, creates a new node otherwise"""
    global unique_table
    
    if not isinstance(x,str):
        if unique_table.__contains__(x):
            return unique_table[x]
        else:
            res=Node(x)
            unique_table[x]=res
        return res
    temp_key=(x,get_int_key(weight1),get_int_key(weight2),node1,node2)
    if unique_table.__contains__(temp_key):
        return unique_table[temp_key]
    else:
        res=Node(x)
        res.out_weight=[weight1,weight2]
        res.successor=[node1,node2]
        unique_table[temp_key]=res
    return res

def find_computed_table(item):
    """To return the results that already exist"""
    global computed_table
    the_keys = []
    if item[0]=='s':
        the_keys.append(('s',get_int_key(item[1].weight),item[1].node,item[2],item[3]))
    elif item[0] == '+':
        the_keys.append(('+',get_int_key(item[1].weight),item[1].node,get_int_key(item[2].weight),item[2].node))
        the_keys.append(('+',get_int_key(item[2].weight),item[2].node,get_int_key(item[1].weight),item[1].node))
    else:
        the_keys.append(('*',get_int_key(item[1].weight),item[1].node,get_int_key(item[2].weight),item[2].node,item[3]))
        the_keys.append(('*',get_int_key(item[2].weight),item[2].node,get_int_key(item[1].weight),item[1].node,item[3]))
    for the_key in the_keys:
        if computed_table.__contains__(the_key):
            res = computed_table[the_key]
            tdd = TDD(res[1])
            tdd.weight = res[0]
            return tdd
    return None

def insert_2_computed_table(item,res):
    """To insert an item to the computed table"""
    global computed_table
    if item[0]=='s':
        the_key = ('s',get_int_key(item[1].weight),item[1].node,item[2],item[3])
    elif item[0] == '+':
        the_key = ('+',get_int_key(item[1].weight),item[1].node,get_int_key(item[2].weight),item[2].node)
    else:
        the_key = ('*',get_int_key(item[1].weight),item[1].node,get_int_key(item[2].weight),item[2].node,item[3])
    computed_table[the_key] = (res.weight,res.node)      

def Matrix2TDD(U,var):
    """Get the TDD of a 2*2 matrix, here var = [column index, row index]"""
    global global_index_order
    (m,n)=U.shape
    if m==n==1:
        node=Find_Or_Add_Unique_table(1,0,0,None,None)
        res=TDD(node)
        res.weight=U[0][0]
        return res

    if m==2 and n==1:
        x=var[0].key
        low=Matrix2TDD(U[0:m//2,:],[])
        high=Matrix2TDD(U[m//2:,:],[])
        tdd = normalize(x, low, high)
        tdd.var_set = copy.copy(var)
        return tdd

    if m==1 and n==2:
        x=var[0].key
        low=Matrix2TDD(U[:,0:n//2],[])
        high=Matrix2TDD(U[:,n//2:],[])
        tdd = normalize(x, low, high)
        tdd.var_set = copy.copy(var)
        return tdd

    if var[0] < var[1]:
        x=var[0].key
        low=Matrix2TDD(U[:,0:n//2],[var[1]])
        high=Matrix2TDD(U[:,n//2:],[var[1]])
    else:
        x=var[1].key
        low=Matrix2TDD(U[0:m//2,:],[var[0]])
        high=Matrix2TDD(U[m//2:,:],[var[0]])
    tdd=normalize(x,low,high)
    tdd.index_set = copy.copy(var)
    return tdd


def diag_matrix_2_TDD(U,var):
    """Get the TDD of a diagonal matrix; var is the index this matrix"""
    node1=Find_Or_Add_Unique_table(1,0,0,None,None)
    low=TDD(node1)
    low.weight=U[0][0]
    node2=Find_Or_Add_Unique_table(1,0,0,None,None)
    high=TDD(node2)
    high.weight=U[1][1]
    if high==low:
        return high
    
    x=var[0].key
    tdd=normalize(x,low,high)
    tdd.index_set = copy.copy(var)
    return tdd

def cnot_2_TDD(var,case=1):
    """Get the TDD of a CNOT gate"""
    global global_index_order
    
    X=np.array([[0,1],[1,0]])
    I=np.array([[1,0],[0,1]])
    
    if case==1:
        if var[0] < var[3] and var[0] < var[4]:
            x=var[0].key
            low=Matrix2TDD(I,[var[3],var[4]])
            high=Matrix2TDD(X,[var[3],var[4]])
            tdd = normalize(x,low,high)
            tdd.index_set = [var[0],var[2],var[3],var[4]]
            return tdd
        elif var[3] < var[0] and var[3] < var[4]:
            x=var[3].key
            low=Matrix2TDD(I,[var[4],var[0]])
            high=Matrix2TDD(X,[var[4],var[0]])
            tdd = normalize(x,low,high)
            tdd.index_set = [var[0],var[2],var[3],var[4]]
            return tdd
        else:
            x = var[4].key
            low = Matrix2TDD(I, [var[3], var[0]])
            high = Matrix2TDD(X, [var[3], var[0]])
            tdd = normalize(x,low,high)
            tdd.index_set = [var[0],var[2],var[3],var[4]]
            return tdd
    if case==2:
        term=Find_Or_Add_Unique_table(1,0,0,None,None)
        tdd=TDD(term)
        tdd.index_set = [var[0],var[1],var[2]]
        return tdd
    
    if case==3:
        if var[1] < var[3] and var[1] < var[4]:
            x=var[1].key
            low=Matrix2TDD(I,[var[3],var[4]])
            high=Matrix2TDD(X,[var[3],var[4]])
            tdd = normalize(x,low,high)
            tdd.index_set = [var[1],var[3],var[4]]
            return tdd
        elif var[3] < var[1] and var[3] < var[4]:
            x=var[3].key
            low=Matrix2TDD(I,[var[4],var[1]])
            high=Matrix2TDD(X,[var[4],var[1]])
            tdd = normalize(x,low,high)
            tdd.index_set = [var[1],var[3],var[4]]
            return tdd
        else:
            x = var[4].key
            low = Matrix2TDD(I, [var[3],var[1]])
            high = Matrix2TDD(X, [var[3],var[1]])
            tdd = normalize(x,low,high)
            tdd.index_set = [var[1],var[3],var[4]]
            return tdd
        
def Slicing(tdd,x,c):
    """Slice a TDD with respect to x=c"""
    global global_index_order
    
    v=tdd.node.key
    
    temp_index_set=[]
    for var in tdd.index_set:
        if var.key !=x:
            temp_index_set.append(var)
            
    temp_index_set2=[]
    for var in tdd.index_set:
        if var.key !=v:
            temp_index_set2.append(var)            
    
    if find_computed_table(['s',tdd,x,c]):
        res = find_computed_table(['s',tdd,x,c])
        res.index_set = temp_index_set
        return res
    
    if not isinstance(v,str) or global_index_order.index(v)>global_index_order.index(x):
        res = copy.copy(tdd)
        res.index_set = temp_index_set
        insert_2_computed_table(['s',tdd,x,c],res)
        return res
    term=Find_Or_Add_Unique_table(1,0,0,None,None)
    leftweight=tdd.node.out_weight[0]*tdd.weight
    if get_int_key(leftweight)==(0,0):
        leftChild=TDD(term)
        leftChild.weight=0
    else:
        leftChild=TDD(tdd.node.successor[0])
        leftChild.weight=leftweight
    rightweight=tdd.node.out_weight[1]*tdd.weight
    if get_int_key(rightweight)==(0,0):
        rightChild=TDD(term)
        rightChild.weight=0
    else:
        rightChild=TDD(tdd.node.successor[1])
        rightChild.weight=rightweight            
    leftChild.index_set=temp_index_set2
    rightChild.index_set=temp_index_set2
    
    if v==x:
        if c==0:
            return leftChild
        else:
            return rightChild

    low=Slicing(leftChild,x,c)
    high=Slicing(rightChild,x,c)
    res=normalize(v,low,high)
    res.index_set = temp_index_set
    insert_2_computed_table(['s',tdd,x,c],res)
    return res

def apply(tdd1,tdd2,op):
    """The apply function of two TDDs. Mostly, it is used to do addition here."""
    global global_index_order    
    
    v1=tdd1.node.key
    v2=tdd2.node.key
    if find_computed_table(['+',tdd1,tdd2]):
        return find_computed_table(['+',tdd1,tdd2]) 
    
    term=Find_Or_Add_Unique_table(1,0,0,None,None)
    if op is '+':
        if tdd1.node==tdd2.node:
            weig=tdd1.weight+tdd2.weight
            if get_int_key(weig)==(0,0):
                res=TDD(term)
                res.weight=0
                return res
            else:
                res=TDD(tdd1.node)
                res.weight=weig
                insert_2_computed_table(['+',tdd1,tdd2],res)
                return res
    if op is '*':
        if v1==1:
            res=TDD(tdd2.node)
            res.weight=tdd1.weight*tdd2.weight
            insert_2_computed_table(['+',tdd1,tdd2],res)
            return res
        if v2==1:
            res=TDD(tdd1.node)
            res.weight=tdd1.weight*tdd2.weight
            insert_2_computed_table(['+',tdd1,tdd2],res)
            return res
    if isinstance(v1,str) and not isinstance(v2,str):
        x=v1        
    elif isinstance(v2,str) and not isinstance(v1,str):
        x=v2    
    elif global_index_order.index(v1)<=global_index_order.index(v2):
        x=v1
    else:
        x=v2
        
    low=apply(Slicing(tdd1,x,0),Slicing(tdd2,x,0),op)
    high=apply(Slicing(tdd1,x,1),Slicing(tdd2,x,1),op)
    res = normalize(x,low,high)
    insert_2_computed_table(['+',tdd1,tdd2],res)
    return res

def contraction(tdd1,tdd2):
    """The contraction of two TDDs"""
    global global_index_order
    v1=tdd1.node.key
    v2=tdd2.node.key
    term=Find_Or_Add_Unique_table(1,0,0,None,None)
    var_out=[]
    var_conted=[]
    for var in tdd1.index_set:
        if var in tdd2.index_set:
            if not var.key in var_conted:
                var_conted.append(var.key)
        else:
            var_out.append(var)
    for var in tdd2.index_set:
        if not var in tdd1.index_set:
            var_out.append(var)
    var_out_key=[]
    for var in var_out:
        var_out_key.append(var.key)
        
    var_temp=[]
    for var in var_conted:
        if not var in var_out_key:
            var_temp.append((global_index_order.index(var),var))
    var_temp.sort()
    var_conted=[]
    for item in var_temp:
        var_conted.append(item[1])
    if find_computed_table(['*',tdd1,tdd2,tuple(var_conted)]):
        tdd = find_computed_table(['*',tdd1,tdd2,tuple(var_conted)])
        tdd.index_set = copy.copy(var_out)
        return  tdd    
    
    if v1==1 and v2==1:
        weig=(2**len(var_conted))*tdd1.weight*tdd2.weight
        tdd=TDD(term)
        if get_int_key(weig)==(0,0):
            tdd.weight=0
        else:
            tdd.weight=weig
        tdd.index_set = copy.copy(var_out)
        
        insert_2_computed_table(['*',tdd1,tdd2,tuple(var_conted)],tdd)
        return tdd

    if v1==1:
        if len(var_conted)==0:
            weig=tdd1.weight*tdd2.weight
            if get_int_key(weig)==(0,0):
                tdd=TDD(term)
                tdd.weight=0
                tdd.index_set = copy.copy(var_out)
                insert_2_computed_table(['*',tdd1,tdd2,tuple(var_conted)],tdd)
                return tdd
            else:
                node=tdd2.node
                tdd=TDD(node)
                tdd.weight=weig
                tdd.index_set = copy.copy(var_out)
                insert_2_computed_table(['*',tdd1,tdd2,tuple(var_conted)],tdd)
                return tdd
        else:
            x = var_conted[0]
            low = contraction(tdd1,Slicing(tdd2,x,0))
            high = contraction(tdd1,Slicing(tdd2,x,1))
            tdd = apply(low,high,'+')
            tdd.index_set = copy.copy(var_out)
            insert_2_computed_table(['*',tdd1,tdd2,tuple(var_conted)],tdd)
            return tdd
            
    if v2==1:
        return contraction(tdd2,tdd1)
        
    if global_index_order.index(v1)<=global_index_order.index(v2):
        x=v1
    else:
        x=v2
        
    if not x in var_conted:
        low=contraction(Slicing(tdd1,x,0),Slicing(tdd2,x,0))
        high=contraction(Slicing(tdd1,x,1),Slicing(tdd2,x,1))
        tdd=normalize(x,low,high)
        tdd.index_set = copy.copy(var_out)
        insert_2_computed_table(['*',tdd1,tdd2,tuple(var_conted)],tdd)
        return tdd
    
    else:
        low=contraction(Slicing(tdd1,x,0),Slicing(tdd2,x,0))
        high=contraction(Slicing(tdd1,x,1),Slicing(tdd2,x,1))
        tdd=apply(low,high,'+')
        tdd.index_set = copy.copy(var_out)
        insert_2_computed_table(['*',tdd1,tdd2,tuple(var_conted)],tdd)
        return tdd
    