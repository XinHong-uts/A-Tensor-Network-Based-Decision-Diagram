from TDD.TDD import Index

def get_tensor_index(cir, num_qubit):
    """return a dict that link every quantum gate to the corresponding index"""
    qubits_index=[0]*num_qubit
    dict1=dict()
    hyper_index=dict()
    for k in cir.nodes():
        operation=cir.nodes[k]['operation']
        nam=operation.name
        if nam != 'CX':
            q=operation.involve_qubits[0]
            var_in='x'+ str(q)+str(0)+str(qubits_index[q])
            var_out='x'+ str(q)+str(0)+str(qubits_index[q]+1)
            if not var_in in hyper_index:
                    hyper_index[var_in]=0
            if not var_out in hyper_index:
                    hyper_index[var_out]=0
            if nam =='z' or nam == 't' or nam == 's'or nam == 'tdg' or nam == 'sdg'or nam == 'rz':
                dict1[k]=[Index(var_in,hyper_index[var_in]),Index(var_in,hyper_index[var_in]+1)]
                hyper_index[var_in]+=1
            else:
                dict1[k]=[Index(var_in,hyper_index[var_in]),Index(var_out,hyper_index[var_out])]
                qubits_index[q]+=1
        else:
            q=[operation.control_qubit,operation.target_qubit]
            var_con='x'+ str(q[0])+str(0)+str(qubits_index[q[0]])
            var_tar_in='x'+ str(q[1])+str(0)+str(qubits_index[q[1]])
            var_tar_out='x'+ str(q[1])+str(0)+str(qubits_index[q[1]]+1)
            if not var_con in hyper_index:
                hyper_index[var_con]=0
            if not var_tar_in in hyper_index:
                hyper_index[var_tar_in]=0
            if not var_tar_out in hyper_index:
                hyper_index[var_tar_out]=0
            cnot_index=[]
            cnot_index.append(Index(var_con,hyper_index[var_con]))
            cnot_index.append(Index(var_con,hyper_index[var_con]+1))
            cnot_index.append(Index(var_con,hyper_index[var_con]+2))
            cnot_index.append(Index(var_tar_in,hyper_index[var_tar_in]))
            cnot_index.append(Index(var_tar_out,hyper_index[var_tar_out]))
            dict1[k]=cnot_index
            hyper_index[var_con]+=2
            qubits_index[q[1]]+=1
        
    for k in range(num_qubit):
        last1='x'+str(k)+str(0)+str(qubits_index[k])
        new1='y'+str(k)
        last2='x'+str(k)+str(0)+str(0)
        new2='x'+str(k)
        for m in dict1:
            dict1[m]=[Index(new1,item.idx) if item.key ==last1 else item for item in dict1[m]]
            dict1[m]=[Index(new2,item.idx) if item.key ==last2 else item for item in dict1[m]]

    return dict1

def get_real_qubit_num(cir):
    qubit_num=0
    for k in cir.nodes():
        temp=max(cir.nodes[k]['operation'].involve_qubits)
        qubit_num=max(qubit_num,temp)
    qubit_num+=1
    return qubit_num

def get_gates_number(cir):
    return len(cir.nodes)

def circuit_partion1(cir,num_qubit,cx_max):
    """The first partition scheme;
    cx_max is the number of CNOT that allowed to be cut every time
    """
    res=[[[],[]]]
    cx_num=0
    level=0
    qubits=[]
        
    qubits = []    
    qubits.append([k for k in range(num_qubit//2)])
    qubits.append([k for k in range(num_qubit//2,num_qubit)])

    for k in cir.nodes():
        opera=cir.nodes[k]['operation']
        nam=opera.name
        if nam!='CX':
            if opera.involve_qubits[0] in qubits[0]:
                res[level][0].append(k)
            else:
                res[level][1].append(k)
        else:
            if opera.control_qubit in qubits[0] and opera.target_qubit in qubits[0]:
                res[level][0].append(k)
            elif opera.control_qubit in qubits[1] and opera.target_qubit in qubits[1]:
                res[level][1].append(k)
            else:
                cx_num+=1
                if cx_num<=cx_max:
                    res[level][0].append(k)
                    res[level][1].append(k)
                else:
                    level+=1
                    res.append([])
                    res[level].append([])
                    res[level].append([])
                    res[level][0].append(k)
                    res[level][1].append(k)
                    cx_num=1
    print('circuit blocks:',2*(level+1))
    return res,qubits

def circuit_partion2(cir,num_qubit,cx_max,c_part_width):
    """The first partition scheme;
    cx_max is the number of CNOT that allowed to be cut every time
    """
    res = [[[],[],[]]]
    cx_num=0
    level=0   
    
    qubits=[]
    qubits.append([k for k in range(num_qubit//2)])
    qubits.append([k for k in range(num_qubit//2,num_qubit)])
    qubits.append([])
    c_range=[num_qubit,0]
    for k in cir.nodes():
        opera=cir.nodes[k]['operation']
        nam=opera.name
        if nam!='CX':
            if opera.involve_qubits[0] in qubits[0]:
                res[level][0].append(k)
            elif opera.involve_qubits[0] in qubits[1]:
                res[level][1].append(k)
            else:
                res[level][2].append(k)
        else:
            if opera.control_qubit in qubits[0] and opera.target_qubit in qubits[0]:
                res[level][0].append(k)
            elif opera.control_qubit in qubits[1] and opera.target_qubit in qubits[1]:
                res[level][1].append(k)
            elif opera.control_qubit in qubits[2] and opera.target_qubit in qubits[2]:
                res[level][2].append(k)
            else:
                if cx_num < cx_max:
                    cx_num+=1
                    res[level][0].append(k)
                    res[level][1].append(k)
                else:
                    c_width=max(c_range[1],max(opera.involve_qubits))-min(c_range[0],min(opera.involve_qubits))+1
                    if c_width<c_part_width:
                        res[level][2].append(k)
                        c_range[0]=min(c_range[0],min(opera.involve_qubits))
                        c_range[1]=max(c_range[1],max(opera.involve_qubits))
                        qubits[0]=[k for k in range(0,c_range[0])]
                        qubits[1]=[k for k in range(c_range[1]+1,num_qubit)]
                        qubits[2]=[k for k in range(c_range[0],c_range[1]+1)]
                    else:
                        level+=1
                        res.append([])
                        res[level].append([])
                        res[level].append([])
                        res[level].append([])
                        qubits.clear()
                        qubits.append([k for k in range(num_qubit//2)])
                        qubits.append([k for k in range(num_qubit//2,num_qubit)])
                        qubits.append([])
                        c_range=[num_qubit,0]
                        if opera.control_qubit in qubits[0] and opera.target_qubit in qubits[0]:
                            res[level][0].append(k)
                            cx_num=0
                        elif opera.control_qubit in qubits[1] and opera.target_qubit in qubits[1]:
                            res[level][1].append(k)
                            cx_num=0
                        else:
                            res[level][0].append(k)
                            res[level][1].append(k)
                            cx_num=1
    print('circuit blocks:',3*(level+1))
    return res
