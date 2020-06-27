from operator import itemgetter, attrgetter
import sys
#sys.path.append('./rins_comp')
from Aminoacid import *

#CRIAÇAO DOS ARRAYS DE OBJETOS
#ENTRADA - ARQUIVOS NODES/EDGES
#SAIDA - VETOR DE OBJETOS(AMINOACIDOS) DE CADA ARQUIVO
def array_construct(directory):
    vet_obj = []
    file = open(directory, 'r')
    #name = directory.replace('/', '_').split('_')[2]
    for line in file.readlines():
        line = line.split("\t")
        if 'Chain' in line or 'Interaction' in line:
            x = line[1]
            continue
        else:
            if x == 'Chain':
                aminoacid_obj = Aminoacid_node(line[0], line[1], int(line[2]), line[3], line[4], line[5], line[6],
                                               line[7], line[8], line[9], line[10], line[11], line[12], line[13])
            else:
                aminoacid_obj = Aminoacid_edge(line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7],
                                               line[8], line[9], line[10], line[11])
            vet_obj.append(aminoacid_obj)
    if x == 'Chain':
        vet_obj = sorted(vet_obj, key=attrgetter('_chain', '_position'))
 #       file_exit = open(directory_exit, 'w')

 #       for i in vet_obj:
 #           file_exit.write(i.show_id())
 #           file_exit.write('\n')
 #       file_exit.close()
        
    else:
        vet_obj = sorted(vet_obj, key=attrgetter('_chain', '_numNode1', '_chain2', '_numNode2'))

    file.close()

    

    return vet_obj


#CRIAÇAO DOS DICIONARIOS
def dict_construct(vet):
    dict = {}
    for i in vet:
        if i._chain in dict:
            dict[i._chain].append(i)
        else:
            dict[i._chain] = [i]
    return dict


#CONSTRUIR DICIONARIOS DE SIMILARIDADE
#ENTRADA - DICIONARIO PROTEIN 1 - DICIONARIO PROTEIN 2
def sim_dict(dict_1, dict_2):
    dict = {}
    for chain_1 in dict_1: #1b55 chain A,B
        for chain_2 in dict_2:  # 2z0p chain A,B,C,D
            soma_dif = 0
            name_key = chain_1+' x '+chain_2
            dict[name_key] = 0
            for i in dict_1[chain_1]: #1b55 node
                ver = True
                for j in dict_2[chain_2]: #2z0p node
                    if i.compare_nodeId(j):
                        ver = False
                        break
                if ver == True:
                    soma_dif = soma_dif + 1
            dict[name_key] = round(((len(dict_1[chain_1]) - soma_dif) * 100)/(len(dict_1[chain_1])), 2)
    return dict


#COMPARAÇAO DOS ARQUIVOS NODES PROTEIN 1 X PROTEIN 2
#ENTRADA - DICIONARIO PROTEIN 1 - DICIONARIO PROTEIN 2
#SAIDA - ARQUIVO DE COMPARAÇAO DE NODES
def compare_nodes(dict_1, dict_2, num, directory):
    dic_sim = sim_dict(dict_1, dict_2)
    file = open(directory, 'w')
    dict_keys_1 = dict_1.keys()
    protein_name_1 = str(dict_keys_1).replace('dict_keys', '').split('\'')[1]
    protein_name_1 = dict_1[protein_name_1][0]._pdbFileName.split('.')[0]
    dict_keys_2 = dict_2.keys()
    protein_name_2 = str(dict_keys_1).replace('dict_keys', '').split('\'')[1]
    protein_name_2 = dict_2[protein_name_2][0]._pdbFileName.split('.')[0]
    for chain_1 in dict_1: #1b55 chain A,B
        for chain_2 in dict_2:  # 2z0p chain A,B,C,D
            key = chain_1 +' x '+ chain_2
            if dic_sim[key] > num:
                #
                file_chain = open(directory + '-' +chain_1 + 'x' + chain_2 + '.txt', 'w')
                #
                soma_dif = 0
                x = protein_name_1 + ' - ' + chain_1 + ' x ' + chain_2 + ' - ' + protein_name_2
                file.write(protein_name_1 + ' - ' + chain_1 + ' x ' + chain_2 + ' - ' + protein_name_2 + '\n')
                file.write('NodeId\t\t\tChain\t\tPosition\tDegree\t\tBfactorCA\t\tpdbFileName\n')
                for i in dict_1[chain_1]: #1b55 node
                    ver = True
                    for j in dict_2[chain_2]: #2z0p node
                        if i.compare_nodeId(j):
                            ver = False
                            break
                    if ver == True:
                        file.write(i.show_amino()+'\n')
                        #
                        file_chain.write(i.show_id_chain()+'\n')
                        #
                        soma_dif = soma_dif + 1
                file.write('Total: '+str(soma_dif)+' Aminoacids \n')
                file.write('Similaridade: ' + str(dic_sim[key]) + '% \n')
                file.write('\n')
                    



#COMPARAÇAO DE ATRIBUTOS - DEGREE OU BFACTOR
#ENTRADA - DICIONARIO PROTEIN 1 - DICIONARIO PROTEIN 2
#SAIDA - ARQUIVO DE COMPARAÇAO DE ATRIBUTO
def compare_atrib(dict_1, dict_2, attrib, directory):
    file = open(directory, 'w')
    dict_keys_1 = dict_1.keys()
    protein_name_1 = str(dict_keys_1).replace('dict_keys', '').split('\'')[1]
    protein_name_1 = dict_1[protein_name_1][0]._pdbFileName.split('.')[0]
    dict_keys_2 = dict_2.keys()
    protein_name_2 = str(dict_keys_1).replace('dict_keys', '').split('\'')[1]
    protein_name_2 = dict_2[protein_name_2][0]._pdbFileName.split('.')[0]
    for chain_1 in dict_1:
        for chain_2 in dict_2:
            soma_dif = 0
            file.write(protein_name_1 + ' - ' + chain_1 + ' x ' + chain_2 + ' - ' + protein_name_2 + '\n')
            file.write('NodeId\t\t\tChain\t\tPosition\tDegree\t\tBfactorCA\t\tpdbFileName\n')
            for i in dict_1[chain_1]:
                for j in dict_2[chain_2]:
                    if attrib == 'degree':
                        if i.compare_nodeId(j) and not i.compare_degree(j):
                            soma_dif = soma_dif + 1
                            file.write(i.show_amino()+ '\n')
                            file.write(j.show_amino()+'\n')
                            #if int(i._degree) > int(j._degree):
                            #    score = int(j._degree) / int(i._degree)
                            #else:
                            #    score = int(i._degree) / int(j._degree)
                            #score = round(score, 2)
                            #file.write('Score: ' + str(score) + '\n')
                            #file.write('\n')
                    elif attrib == 'bfactor':
                        if i.compare_nodeId(j) and not i.compare_bfactor(j):
                            soma_dif = soma_dif + 1
                            file.write(i.show_amino()+ '\n')
                            file.write(j.show_amino()+'\n')
            file.write('Total: ' + str(soma_dif) + ' Aminoacids\n')
            file.write('\n')


#COMPARAÇAO DOS ARQUIVOS EDGES PROTEIN 1 X PROTEIN 2
#ENTRADA - DICIONARIO PROTEIN 1 - DICIONARIO PROTEIN 2
#SAIDA - ARQUIVO DE COMPARAÇAO DE EDGES
def compare_edges(dict_1, dict_2, num, directory):
    vet_aux1 = []
    vet_aux2 = []
    file = open(directory, 'w')
  
    for chain_1 in dict_1: #1b55 chain A,B
        vet_aux1 = []

        for x in dict_1[chain_1]:
            vet_aux1.append(x)

        for chain_2 in dict_2:  # 2z0p chain A,B,C,D
            vet_aux2 = []

            for y in dict_2[chain_2]:
                vet_aux2.append(y)

            file.write('NodeId1 - ' + chain_1+'\t\t\t Interaction \t\t '+chain_2+' - NodeId2\n')

            for i in vet_aux1: #1b55 node
                for j in vet_aux2: #2z0p node                    
                    if i.compare_edge(j):                      
                        vet_aux2.remove(j)
                        break

            for w in vet_aux2:
                file.write(w.show_edge() + '\n')

