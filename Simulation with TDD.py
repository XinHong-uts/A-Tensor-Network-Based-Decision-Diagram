import time
from TDD.TDD import Ini_TDD,TDD,Find_Or_Add_Unique_table
from TDD.TDD import Matrix2TDD,diag_matrix_2_TDD,cnot_2_TDD,contraction
# from TDD.TDD_show import TDD_show # This is just used to Graphical display the qmdd
from cir_input.qasm import CreateCircuitFromQASM
from cir_input.circuit_DG import CreateDGfromQASMfile
from cir_input.circuit_process import get_tensor_index,get_real_qubit_num,get_gates_number
from cir_input.circuit_process import circuit_partion1,circuit_partion2
from func_timeout import func_set_timeout
import os

def get_tdd(operation,var_list,involve_qubits):
    """get the TDD of the correct part of quantum gate"""
    nam=operation.name
    if nam =='z' or nam == 't' or nam == 's'or nam == 'tdg' or nam == 'sdg'or nam == 'rz':
        return diag_matrix_2_TDD(operation.u_matrix,var_list)
    elif nam =='CX':
        if operation.control_qubit in involve_qubits and operation.target_qubit in involve_qubits:
            return cnot_2_TDD(var_list,case=1)
        if operation.control_qubit in involve_qubits and not operation.target_qubit in involve_qubits:
            return cnot_2_TDD(var_list,case=2)
        else:
            return cnot_2_TDD(var_list,case=3)
    else:
        return Matrix2TDD(operation.u_matrix,var_list)
    
def get_tdd_of_a_part_circuit(involve_nodes,involve_qubits,cir,node_var):
    """get the TDD of a part of circuit"""
    node=Find_Or_Add_Unique_table(1,0,0,None,None)
    tdd=TDD(node)
    max_node_num = 0
    max_node_num = max(max_node_num,tdd.node_number())
    for k in involve_nodes:
        temp_tdd=get_tdd(cir.nodes[k]['operation'],node_var[k],involve_qubits)
        tdd=contraction(tdd,temp_tdd)
        max_node_num = max(tdd.node_number(),max_node_num)
        
    return tdd,max_node_num

@func_set_timeout(3600)
def Simulation_with_cir_partion(cir, num_qubit,partion_scheme=0,cx_max=2,c_part_width=2):
    """Simulate a circuit with TDD;
    cir is a dag circuit;
    num_qubit is the real qubit number of the circuit;
    partition_scheme = 0, 1, 2
    cx_max is the number of CNOT that allowed to be cut in two partition scheme
    c_part_width is the allowed qubit numbers of the C part in the second partition scheme;
    return the result tdd and other corresponding information
    """
    var_order=[]
    gate_num=len(cir.nodes)
    for k in range(num_qubit-1,-1,-1):
        var_order.append('x'+str(k))
        for j in range(gate_num):
            var_order.append('x'+str(k)+str(0)+str(j))
        var_order.append('y'+str(k))
    
    node_var = get_tensor_index(cir, num_qubit)
    
    max_node_num = 0

    Ini_TDD(var_order)
    if partion_scheme==0:
        involve_qubits=[k for k in range(num_qubit)]
        tdd,max_node_num = get_tdd_of_a_part_circuit(cir.nodes,involve_qubits,cir,node_var)
        return tdd, max_node_num,1
    
    
    if partion_scheme==1:
        partion_cir,involve_qubits=circuit_partion1(cir,num_qubit,cx_max)
        node = Find_Or_Add_Unique_table(1,0,0,None,None)
        tdd = TDD(node)
        max_node_num = max(max_node_num,tdd.node_number())
        for level in range(len(partion_cir)):
            tdd1,node_num1=get_tdd_of_a_part_circuit(partion_cir[level][0],involve_qubits[0],cir,node_var)
            max_node_num = max(max_node_num,node_num1)
            
            tdd2,node_num2=get_tdd_of_a_part_circuit(partion_cir[level][1],involve_qubits[1],cir,node_var)
            max_node_num = max(max_node_num,node_num2)
                                   
            temp_tdd=contraction(tdd1,tdd2)
            max_node_num = max(max_node_num,temp_tdd.node_number())

            tdd=contraction(tdd,temp_tdd)                                 
            max_node_num = max(max_node_num,tdd.node_number())

        return tdd, max_node_num,2*len(partion_cir)
    
    if partion_scheme==2:
        involve_qubits=[]
        involve_qubits.append([k for k in range(num_qubit//2)])
        involve_qubits.append([k for k in range(num_qubit//2,num_qubit)])
        involve_qubits.append([k for k in range(num_qubit)])
        partion_cir=circuit_partion2(cir,num_qubit,cx_max,c_part_width)
        node=Find_Or_Add_Unique_table(1,0,0,None,None)
        tdd=TDD(node)
        max_node_num = max(max_node_num,tdd.node_number())
        for level in range(len(partion_cir)):
            tdd1,node_num1=get_tdd_of_a_part_circuit(partion_cir[level][0],involve_qubits[0],cir,node_var)
            max_node_num = max(max_node_num,node_num1)

            tdd2,node_num2=get_tdd_of_a_part_circuit(partion_cir[level][1],involve_qubits[1],cir,node_var)
            max_node_num = max(max_node_num,node_num2)

            temp_tdd=contraction(tdd1,tdd2)
            max_node_num = max(max_node_num,temp_tdd.node_number())
            
            tdd3,node_num3=get_tdd_of_a_part_circuit(partion_cir[level][2],involve_qubits[2],cir,node_var)
            max_node_num = max(max_node_num,node_num3)
            
            temp_tdd=contraction(tdd3,temp_tdd)
            max_node_num = max(max_node_num,temp_tdd.node_number())
            
            tdd=contraction(tdd,temp_tdd)
            max_node_num = max(max_node_num,tdd.node_number())
            
        return tdd,max_node_num,3*len(partion_cir)

    
def TDD_simulate_test(partion_scheme=0):
    path='Benchmarks/'
    file_list = os.listdir(path)
    for file_name in file_list:
        print('circuit：',file_name)
        try:
            cir,res = CreateDGfromQASMfile(file_name, path, flag_single=True)
        except:
            continue
        dag_cir=res[0]

        num_qubit = get_real_qubit_num(dag_cir)
        print('qubits:',num_qubit)
        gate_num = get_gates_number(dag_cir)
        print('gates number:',gate_num)
   
        try:
            t_start = time.time()
            tdd,max_node_num,block_num=Simulation_with_cir_partion(dag_cir, num_qubit,partion_scheme,num_qubit//2,num_qubit//2+1)
            run_time=time.time()-t_start
            print('Run time:',run_time)
            print('Max node num:',max_node_num)
            print('Final node number:',tdd.node_number())
        except:
            print('Time out!')
        print('----------------')
    print('=================')

if __name__=="__main__":
    TDD_simulate_test(partion_scheme=0)
    #partion_scheme=0, is to simulate without partition
    TDD_simulate_test(partion_scheme=1)
    #partion_scheme=1，is corresponding to our first partition scheme
    TDD_simulate_test(partion_scheme=2)
    #partion_scheme=2，is corresponding to our second partition scheme