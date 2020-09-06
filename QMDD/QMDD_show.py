"""This page is just to get the Graphical output of QMDD"""
from graphviz import Digraph, nohtml

from IPython.display import Image

def QMDD_show(qmdd):
    get_idx(qmdd)
    edge=[]              
    dot=Digraph(name='reduced_tree')
    dot=layout(qmdd,dot,edge)
    dot.node('0','',shape='none')
    dot.edge('0',qmdd.node.idx,color="blue",label=str(round(qmdd.weight,2)))
    dot.format = 'png'
    return Image(dot.render('output'))


        
def get_idx(qmdd,idx=1):
    qmdd.node.idx=str(idx)
    idx+=1
    for k in range(4):
        if qmdd.node.successor[k]:
            idx=get_idx(qmdd.node.successor[k],idx)
    return idx
    
def layout(qmdd,dot=Digraph(),succ=[]):
    dot.node(qmdd.node.idx, str(qmdd.node.key), fontname="helvetica",shape="circle",color="red")
        
    col=['red','blue','black','green']
    for k in range(4):
        if qmdd.node.successor[k]:
            if not qmdd.node.successor[k].node.idx in succ:
                layout(qmdd.node.successor[k],dot,succ)
                dot.edge(qmdd.node.idx,qmdd.node.successor[k].node.idx,color=col[k],label=str(round(qmdd.node.successor[k].weight,2)))
            else:
                dot.edge(qmdd.node.idx,qmdd.node.successor[k].node.idx,color=col[k],label=str(round(qmdd.node.successor[k].weight,2)))
            succ.append(qmdd.node.successor[k].node.idx)
    return dot
    
        