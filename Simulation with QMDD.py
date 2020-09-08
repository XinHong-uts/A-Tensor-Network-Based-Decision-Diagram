import time
import numpy as np
from QMDD.QMDD import Ini_QMDD, QMDD_Mul,get_cnot_QMDD_n,get_single_QMDD
# from QMDD.QMDD_show import QMDD_show # This is just used to Graphical display the qmdd
from qiskit.quantum_info.operators import Operator
from cir_input.qasm import CreateCircuitFromQASM
from cir_input.circuit_DG import CreateDGfromQASMfile
from cir_input.circuit_process import get_real_qubit_num,get_gates_number
from func_timeout import func_set_timeout
import os

def CirData2QMDD(data, H=False):
    """get the qmdd from the circuit data"""
    qiskit_gate = data[0]
    qubits = []
    for q in data[1]:
        qubits.append(q.index)
    if qiskit_gate.name != 'cx':
        u_matrix = Operator(qiskit_gate).data
        if H == True: u_matrix = u_matrix.conj().T
        qmdd = get_single_QMDD(u_matrix, qubits)
    else:
        qmdd = get_cnot_QMDD_n(qubits)       
    return qmdd

@func_set_timeout(3600)
def Simulation_with_QMDD(cir,num_qubit):
    """To simulate a quantum circuit with QMDD"""
    data_list1 = cir._data
    ini_qmdd = Ini_QMDD(num_qubit)
    max_node_num = 0
    qmdd_current=get_single_QMDD(np.eye(2), 0)
    
    max_node_num = max(qmdd_current.node_number(),max_node_num)
    while len(data_list1) > 0 :
            qmdd1=CirData2QMDD(data_list1.pop(0))
            max_node_num = max(max_node_num,qmdd1.node_number())
            qmdd_current = QMDD_Mul(qmdd1, qmdd_current)
            max_node_num = max(max_node_num,qmdd_current.node_number())

    return qmdd_current,max_node_num

def QMDD_simulate_test():
    path='Benchmarks/'

    file_list = os.listdir(path)
    for file_name in file_list:
        print('circuit:',file_name)
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
            qmdd, max_node_num=Simulation_with_QMDD(cir, num_qubit)
            run_time=time.time()-t_start
            print('Run time:',run_time)
            print('Max node num:',max_node_num)
            print('Final node number:',qmdd.node_number())
        except:
            print('Time out!')
        print('----------------')
    

if __name__=="__main__":
    QMDD_simulate_test()
