#!/usr/bin/env python
# coding: utf-8




import pandas as pd
import numpy as np
import os
import re
import sys
import subprocess
from subprocess import PIPE, Popen
from zipfile import ZipFile
from math import sqrt
import matplotlib.pyplot as plt


pd.options.mode.chained_assignment = None

def script(directory, project_path):

    print ("ddddddddd", directory) # queria saber pq isso n ta sendo imprimido

    #ESPECIFICAR DIRETORIOS DE ENTRADA E SAIDA
    #  entrada e saida recebe a mesma coisa
    input_directory= directory
    output_directory= directory 

    print('eeeee',  os.listdir(input_directory))
    #DESCOMPACTANDO E REMOVENDO OS .ZIPs
    for file in os.listdir(input_directory):
    #    print(file)
        if file.endswith('.zip'):
            path = input_directory + '/' +file
            file = ZipFile(path, 'r')
            file.extractall(input_directory)
            file.close()
            os.remove(path)

    arq = os.listdir(directory)

    #VERIFICAÇÃO E SUBSTITUIÇAO DO CARACTERE '-' NO NOME DOS ARQUIVOS
    for i in arq:
        aux = i.replace('-', '_')
        os.rename(input_directory+i, input_directory+aux)
    param = os.listdir(input_directory) 



    #CRIANDO DICIONARIO DE ARQUIVOS DE NODES E EDGES
    dict_protein_nodes_files = {}
    dict_protein_edges_files = {}


    def populate_chain(chain):
        aux = str(chain).split(":")
        chain = aux[0]
        return chain

    def populate_aux_target(aux_chain):
        aux = str(aux_chain).split(":")
        aux_target = ':' + aux[1] + ':' + aux[2] + ':' + aux[3]
        return aux_target

    def populate_residue(NodeId):
        aux = str(NodeId).split(":")
        residue = aux[3]
        return residue

    arq_name = []
    for i in param:
        
        arq_name.clear()
        arq_name.append(i[:-10])
        arq_name.append(i[-9:-4])
        
        if arq_name[1] == "nodes":
            
            #LER ARQUIVO E CRIAR DATAFRAME
            dict_protein_nodes_files[arq_name[0]] = pd.read_csv(input_directory+arq_name[0]+ "_nodes.txt", sep="\t", header= 0)
            
            #CRIAR NOVA COLUNA AUX_NODEID COM NODEID
            if 'Accessibility' not in dict_protein_nodes_files[arq_name[0]].columns:
                dict_protein_nodes_files[arq_name[0]].insert(13, 'Accessibility', 0)
            dict_protein_nodes_files[arq_name[0]]['aux_nodeId'] = '-'
            dict_protein_nodes_files[arq_name[0]]['aux_degree'] = '-'
            
            #POPULAR A COLUNA AUX_NODEID
            dict_protein_nodes_files[arq_name[0]]['aux_nodeId'] = dict_protein_nodes_files[arq_name[0]]['NodeId'].apply(populate_aux_target)
            
            #POPULAR A COLUNA AUX_DEGREE
            dict_protein_nodes_files[arq_name[0]]['aux_degree'] = dict_protein_nodes_files[arq_name[0]]['NodeId'].apply(populate_aux_target)
            dict_protein_nodes_files[arq_name[0]]['aux_degree'] = dict_protein_nodes_files[arq_name[0]].apply(lambda row: row.aux_degree + "_D" + str(row.Degree), axis=1)
            
            #POPULAR A COLUNA RESIDUE (ENCONTRADO CASOS DE PDBs COM RESIDUE NAN)
            dict_protein_nodes_files[arq_name[0]]['Residue'] = dict_protein_nodes_files[arq_name[0]]['NodeId'].apply(populate_residue)

        if arq_name[1] == "edges":
            #LER ARQUIVO E CRIAR DATAFRAME
            dict_protein_edges_files[arq_name[0]] = pd.read_csv(input_directory+arq_name[0]+ "_edges.txt", sep="\t", header= 0)
            
            #POPULAR AS COLUNAS AUX_CHAIN 1, 2
            dict_protein_edges_files[arq_name[0]]['aux_chain_1'] = dict_protein_edges_files[arq_name[0]]['NodeId1'].apply(populate_chain)
            dict_protein_edges_files[arq_name[0]]['aux_chain_2'] = dict_protein_edges_files[arq_name[0]]['NodeId2'].apply(populate_chain)



    #CRIANDO DICIONARIOS DE CHAINS
    dict_protein_nodes_chain = {}

    for i in dict_protein_nodes_files:
        chains = pd.unique(dict_protein_nodes_files[i]["Chain"])
        df = dict_protein_nodes_files[i]
        quant_chain= pd.unique(df["Chain"])
        
        for j in quant_chain:
            key_dict = i + '('+ j + ')'
            dict_protein_nodes_chain[key_dict] = (df[df.Chain == j])
        
        

    dict_protein_edges_chain = {}

    for i in dict_protein_edges_files:
        chains = pd.unique(dict_protein_edges_files[i]["aux_chain_1"])
        df = dict_protein_edges_files[i]
        quant_chain= pd.unique(df["aux_chain_1"])

        for j in quant_chain:
            key_dict = i + '('+ j + ')'
            dict_protein_edges_chain[key_dict] = (df[df.aux_chain_1 == j])



    #NODES

    #NODES PRESENTES APENAS NA PROTEINA 1
    dict_result_nodes_1 = {}

    #NODES PRESENTES APENAS NA PROTEINA 2
    dict_result_nodes_2 = {}

    #POSIÇOES QUE SOFRERAM MUDANÇAS DE AMINOACIDOS
    dict_result_nodes_change = {}

    #AMINOACIDOS QUE SOFRERAM MUDANÇAS DE DEGREE - (INCLUI O DEGREE DE AMINO_CHANGE)
    dict_result_degree_change = {}

    #AMINOACIDOS PRESENTES EM AMBAS AS PROTEINAS
    dict_equal_nodes = {}

    vet = []
    for i in dict_protein_nodes_chain:
        for j in dict_protein_nodes_chain:
            if(i != j) and (j + '-' + i not in vet):
                
                dif = dict_protein_nodes_chain[i].merge(dict_protein_nodes_chain[j], on="aux_nodeId", how='outer', suffixes=['', '_'], indicator=True)       
                dif_nodes_protein_1 = dif[dif["_merge"] == "left_only"]
                
                dif_nodes_protein_1 = dif_nodes_protein_1.drop(dif.columns[[14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31]], axis=1)
                dict_result_nodes_1[i+'-'+j] = dif_nodes_protein_1 
                
                equal_nodes = dif[dif["_merge"] == "both"]
                dict_equal_nodes[i+'-'+j] = equal_nodes
                
                dif_nodes_protein_2 = dif[dif["_merge"] == "right_only"]
                dif_nodes_protein_2 = dif_nodes_protein_2.drop(dif.columns[[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]], axis=1)
                dif_nodes_protein_2.columns = dict_protein_nodes_chain[j].columns
                dif_nodes_protein_2.drop(columns=['aux_nodeId', 'aux_degree'], inplace = True)
                dict_result_nodes_2[i+'-'+j] = dif_nodes_protein_2 
                
                
                #POSIÇOES QUE SOFRERAM MUDANÇAS DE AMINOACIDOS     
                amino_change = pd.DataFrame(columns=dif_nodes_protein_1.columns)
                dict_result_nodes_change[i+"-"+j] = amino_change 
                for w, row in dif_nodes_protein_1.iterrows():
                    if row["Position"] in dif_nodes_protein_2["Position"].values:
                        amino_change = amino_change.append(row)
                        aux = dif_nodes_protein_2[dif_nodes_protein_2.Position == row["Position"]]
                        amino_change = amino_change.append(aux)
                        dict_result_nodes_change[i+"-"+j] = amino_change     
                        
                
                #AMINOACIDOS QUE SOFRERAM MUDANÇAS DE DEGREE
                dif_degree = dict_protein_nodes_chain[i].merge(dict_protein_nodes_chain[j], on="aux_degree", how='outer', suffixes=['', '_'], indicator=True)
                
                dif_degree_protein_1 = dif_degree[dif_degree["_merge"] == "left_only"]
                dif_degree_protein_1 = dif_degree_protein_1.drop(dif_degree.columns[[14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31]], axis=1)
                
                dif_degree_protein_2 = dif_degree[dif_degree["_merge"] == "right_only"]
                dif_degree_protein_2 = dif_degree_protein_2.drop(dif_degree.columns[[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]], axis=1)
                dif_degree_protein_2.columns = dict_protein_nodes_chain[j].columns
                dif_degree_protein_2.drop(columns=['aux_nodeId', 'aux_degree'], inplace = True)
                
                degree_change = pd.DataFrame(columns=dif_degree_protein_1.columns)
                dict_result_degree_change[i+"-"+j] = degree_change 
                for w, row in dif_degree_protein_1.iterrows():
                    if row["Position"] in dif_degree_protein_2["Position"].values:
                        degree_change = degree_change.append(row)
                        aux = dif_degree_protein_2[dif_degree_protein_2.Position == row["Position"]]
                        degree_change = degree_change.append(aux)
                        dict_result_degree_change[i+"-"+j] = degree_change 
                vet.append(i + '-' + j)
                
    vet.clear()            



    #EDGES

    #CRIAR DICIONARIO DE CONTAGEM DE EDGES
    dict_cont_edges = {}

    def populate_aux_position(node_id1):
        aux = str(node_id1).split(':')
        return(int(aux[1]))

    def populate_aux_target(aux_chain):
        aux = str(aux_chain).split(":")
        aux_target = ':' + aux[1] + ':' + aux[2] + ':' + aux[3]
        return aux_target

    def populate_interaction(interaction):
        return(interaction)


    for i in dict_protein_edges_chain:

        
        #CONTAR INTERAÇÕES REPETIDAS
        aux = dict_protein_edges_chain[i].pivot_table(index=['NodeId1', 'Interaction', 'NodeId2', 'aux_chain_1', 'aux_chain_2'], aggfunc='size')
        test=pd.DataFrame(aux)
        check=test.reset_index(inplace=True)
        dict_cont_edges[i] = test
        
        #RENOMEAR COLUNAS
        dict_cont_edges[i].columns = ['NodeId1', 'Interaction', 'NodeId2', 'aux_chain_1', 'aux_chain_2', 'Quant']
        
        #CRIAR AUX_NODEID_1 E 2, INTERACTION
        dict_cont_edges[i]['aux_nodeId_1'] = dict_cont_edges[i]['NodeId1'].apply(populate_aux_target)
        dict_cont_edges[i]['aux_nodeId_2'] = dict_cont_edges[i]['NodeId2'].apply(populate_aux_target)
        dict_cont_edges[i]['aux_interaction'] = dict_cont_edges[i]['Interaction'].apply(populate_interaction)
        #ORDENAR APENAS PELO NODEID1
        dict_cont_edges[i]['aux_position'] = dict_cont_edges[i]['NodeId1'].apply(populate_aux_position)
        dict_cont_edges[i] = dict_cont_edges[i].sort_values('aux_position')



    #OS RESULTADOS VÃO ESTAR EM RAZÃO DA QUANTIDADE. CADA INTERAÇÃO TEM O CAMPO QUANTIDADE, 
    #QUE INDICA A QUANTIDADE DAQUELA INTERAÇÃO

    #EDGES PRESENTES APENAS NA PROTEINA 1
    dict_result_edges_1 = {}

    #EDGES PRESENTES APENAS NA PROTEINA 2
    dict_result_edges_2 = {}

    #MERGE DOS DICIONARIOS DE CONTAGEM DE EDGES
    for i in dict_cont_edges:
        for j in dict_cont_edges:
            if(i != j) and (j + '-' + i not in vet):
                dif = dict_cont_edges[i].merge(dict_cont_edges[j], on=['aux_nodeId_1','aux_interaction','aux_nodeId_2'], how='outer', suffixes=['', '_'], indicator=True)
                
                #APENAS PRESENTE NA PROTEINA 1 - 1B55.A
                dif_edges_protein_1 = dif[dif["_merge"] == "left_only"]
                dif_edges_protein_1 = dif_edges_protein_1.drop(dif.columns[[6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17]], axis=1)
                
                #APENAS PRESENTE NA PROTEINA 2 - 2Z0P.A
                dif_edges_protein_2 = dif[dif["_merge"] == "right_only"]
                dif_edges_protein_2 = dif_edges_protein_2.drop(dif.columns[[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 17]], axis=1)
                
                #PRESENTE EM AMBAS AS PROTEINAS
                dif_edges_protein_3 = dif[dif["_merge"] == "both"]
                
                if not dif_edges_protein_3.empty:

                    #REALOCAR AS INTERAÇÕES PARA CORRIGIR A QUANTIDADE
                    dif_edges_protein_3['dif_quant'] = dif_edges_protein_3.apply(lambda x: x['Quant'] - x['Quant_'], axis=1)

                    maiorquezero = dif_edges_protein_3[dif_edges_protein_3.dif_quant > 0]
                    menorquezero = dif_edges_protein_3[dif_edges_protein_3.dif_quant < 0]


                    maiorquezero.drop(columns=['NodeId1_', 'Interaction_','NodeId2_', 'aux_position_', '_merge', 'aux_nodeId_1','aux_nodeId_2', 'aux_chain_1_','aux_chain_2_','Quant_', 'aux_interaction'], inplace = True) 
                    menorquezero.drop(columns=['NodeId1', 'Interaction','NodeId2','aux_chain_1', 'aux_chain_2', 'Quant', 'aux_nodeId_1', 'aux_nodeId_2', 'Quant', 'aux_position', '_merge','aux_interaction'], inplace = True) 



                    if not maiorquezero.empty:
                        #INSERIR A ROW NO DATAFRAME CORRETO DE ACORDO COM A QUANTIDADE - PROTEINA 1
                        for w, row in maiorquezero.iterrows():
                            row.Quant = int(row.dif_quant)
                            dif_edges_protein_1 = dif_edges_protein_1.append(row)

                    
                    if not menorquezero.empty:               
                        #INSERIR A ROW NO DATAFRAME CORRETO DE ACORDO COM A QUANTIDADE - PROTEINA 2
                        for w, row in menorquezero.iterrows():
                            row.Quant_ = abs(int(row.dif_quant))
                            dif_edges_protein_2 = dif_edges_protein_2.append(row)

                            
                dif_edges_protein_1 = dif_edges_protein_1.sort_values('aux_position')
                
                dif_edges_protein_2 = dif_edges_protein_2.sort_values('aux_position_')
                

                if 'dif_quant'in dif_edges_protein_2.columns:
                    dif_edges_protein_2.drop(columns=['dif_quant'], inplace = True)
                    dif_edges_protein_2.columns = ['NodeId1', 'Interaction', 'NodeId2','aux_chain_1', 'aux_chain_2', 'Quant', 'aux_position']
                else:
                    dif_edges_protein_2.columns = ['NodeId1', 'Interaction', 'NodeId2','aux_chain_1', 'aux_chain_2', 'Quant', 'aux_position']

                dict_result_edges_1[i+'-'+j] = dif_edges_protein_1 
                dict_result_edges_2[i+'-'+j] = dif_edges_protein_2 
                
                vet.append(i + '-' + j)
                            

    vet.clear()      



    #SAIDA DE ARQUIVOS

    #REPORT NODES
    for i in dict_result_nodes_1.keys():
        os.mkdir(output_directory+'/'+i)
        file = open(output_directory+"/"+i+"/report_dif_nodes_"+i+".txt", "a", newline='')
        file.write(i + "\n\n")
        protein_1_name = i.split("-")[0]
        protein_2_name = i.split("-")[1]
        file.write("- AMINOACIDS ONLY PRESENT IN PROTEIN -> " + protein_1_name + "\n\n")
        dict_result_nodes_1[i].to_csv(file, sep="\t", index=False)
        file.write("\n- AMINOACIDS ONLY PRESENT IN PROTEIN -> " + protein_2_name + "\n\n")
        dict_result_nodes_2[i].to_csv(file, sep="\t", index=False)
        file.write("\n- POSITIONS THAT HAVE CHANGED AMINOACIDS -> " + i + "\n\n")
        dict_result_nodes_change[i].to_csv(file, sep="\t", index=False)
        file.close()
        
    #REPORT NODES    
    for i in dict_result_degree_change.keys():
        file = open(output_directory+"/"+i+"/report_dif_degree_"+i+".txt", "a", newline='')
        dict_result_degree_change[i].to_csv(file, sep="\t", index=False)
        file.close()
        
    #REPORT EDGES
    for i in dict_result_edges_1.keys():
        file = open(output_directory+"/"+i+"/report_dif_edges_"+i+".txt", "a", newline='')
        file.write(i + "\n\n")
        protein_1_name = i.split("-")[0]
        protein_2_name = i.split("-")[1]
        file.write("- EDGES ONLY PRESENT IN PROTEIN -> " + protein_1_name + "\n\n")
        dict_result_edges_1[i].to_csv(file, sep="\t", index=False)
        file.write("\n- EDGES ONLY PRESENT IN PROTEIN -> " + protein_2_name + "\n\n")
        dict_result_edges_2[i].to_csv(file, sep="\t", index=False)
        file.close()



    #GERAÇAO DE ARQUIVOS DE GRÁFICOS

    #ORDENAR DATAFRAMES
    def populate_position(position):
        return(int(position))

    def populate_aux_ord(node_id1):
        aux = str(node_id1).split(':')
        return(int(aux[1]))

    for i in dict_protein_nodes_chain:
        dict_protein_nodes_chain[i]['Position'] = dict_protein_nodes_chain[i]['Position'].apply(populate_position)
        dict_protein_nodes_chain[i] = dict_protein_nodes_chain[i].sort_values('Position')
        
        dict_protein_edges_chain[i]['aux_ord'] = dict_protein_edges_chain[i]['NodeId1'].apply(populate_aux_ord)
        dict_protein_edges_chain[i] = dict_protein_edges_chain[i].sort_values('aux_ord')
        
        
        
    #GRAFICO HEATMAP DEGREE
    for i in dict_protein_nodes_chain:
        heatmap_degree = dict_protein_nodes_chain[i]
        
        file = open(output_directory+"/heatmap_degree_"+i+".csv", "a", newline='')
        file.write("group,variable,idnode,degree,color\n") 
        
        maior = heatmap_degree['Position'].max()
        
        #DATAFRAME AUXILIAR PARA VERIFICAR ESPAÇOS EM BRANCO
        cont_df_aux = pd.DataFrame()
        cont_df_aux['Position'] = ''
        cont = 1
        for w in range(0, int(maior)):
            cont_df_aux.loc[w, 'Position'] = cont
            cont+=1
        ##
        
        heatmap_degree = heatmap_degree.merge(cont_df_aux, on="Position", how='outer', suffixes=['', '_'], indicator=True)
            
        for w, row in heatmap_degree.iterrows():
            if row._merge == 'right_only':
                heatmap_degree.loc[w,'Degree'] = 0  
                heatmap_degree.loc[w,'NodeId'] = '-' 
        heatmap_degree = heatmap_degree.sort_values('Position')  
        col = 0
        lin = 1
        #dim_max = 9
        dim_max = int(sqrt(maior))
        pos = 1
        cont = 0
        
        for w, row in heatmap_degree.iterrows():
            cont+=1
            degree = 0
            color = int(row['Degree']) * 1.9 + 75


            file.write(str(col) + ',' + str(lin) + ',' + str(row['NodeId']) + ',' + str(row['Degree']) + ',' + str(color) + '\n')
            pos+=1
            if col <= dim_max:
                col+=1
            else:
                lin+=1
                col = 0
            if col == dim_max:
                lin = cont + 1
                col = 0
        file.close()

                
    #GRAFICO DIF_DEGREE
    for i in dict_equal_nodes:
        #VERIFICA SE EXISTE NODE IGUAL
        if not dict_equal_nodes[i].empty:
            dict_equal_nodes[i]['dif_degree'] = dict_equal_nodes[i].apply(lambda x: abs(x['Degree'] - x['Degree_']), axis=1)
            dif_degree_ord = dict_equal_nodes[i].sort_values('dif_degree', ascending=[False]).head(10)
        
        file = open(output_directory+"/"+i+"/graphic_dif_degree_"+i+".csv", "a", newline='')
        file.write('nodes\tdifferences\n')

        #VERIFICA SE EXISTE NODE IGUAL
        if not dict_equal_nodes[i].empty:        
            for w, row in dif_degree_ord.iterrows():
                file.write(str(row['NodeId']) + '\t' + str(row['dif_degree']) + '\n')
        file.close()   
            



    #GRAFICO DIF_INTERACTION 
    #Gráfico que mostra a diferença na quantidade de interações da comparação. 
    #Apenas interações exclusivas da primeira subtraido das interações exclusivas da segunda

    for i in dict_result_edges_1:
        #TOTAL DE INTERAÇÕES PRESENTES APENAS NA 1B55.A
        aux_1 = dict_result_edges_1[i]
        for w, row in aux_1.iterrows():
            if row.Quant > 1:
                for j in range(1, int(row.Quant)):
                    aux_1 = aux_1.append(row, ignore_index=True)

        cont_int_1_series = aux_1.Interaction.value_counts()

        #TOTAL DE INTERAÇÕES PRESENTES APENAS NA 2Z0P.A
        aux_2 = dict_result_edges_2[i]
        for w, row in aux_2.iterrows():
            if row.Quant > 1:
                for j in range(1, int(row.Quant)):
                    aux_2 = aux_2.append(row, ignore_index=True)

        cont_int_2_series = aux_2.Interaction.value_counts()

        #UNIAO DOS NOMES DAS COLUNAS
        index_list = list(set().union(cont_int_1_series.index, cont_int_2_series.index)) 
        cont_int = pd.DataFrame()
        cont_int = cont_int.assign(**dict.fromkeys(index_list, 0))

        cont_int = cont_int.append(cont_int_1_series)
        cont_int = cont_int.append(cont_int_2_series)
        cont_int.reset_index(drop=True, inplace=True)
        cont_int = cont_int.fillna(0)
        cont_int = abs(cont_int.diff(-1))
        cont_int = cont_int.drop([1])
        cont_int = cont_int.sort_values(by=0, axis=1, ascending=[False])
        cont_int = cont_int.transpose()

        #cont_int['id2'] = cont_int.index
        #cont_int.set_index(0,inplace=True)
        #cont_int.T

        file = open(output_directory+"/"+i+"/graphic_dif_interaction_"+i+".csv", "a", newline='')
        file.write("interactions\tlosses\n")
        cont_int.to_csv(file, sep="\t")
        file.close()


        #REMOVER O 0
        file = open(output_directory+"/"+i+"/graphic_dif_interaction_"+i+".csv", "r")
        contents = file.readlines()
        file.close()

        contents.pop(1) # remove the line item from list, by line number, starts from 0

        file = open(output_directory+"/"+i+"/graphic_dif_interaction_"+i+".csv", "w")
        contents = "".join(contents)
        file.write(contents)
        file.close()


    #GRAFICO HEATMAP - número de diferenças dos edges para cada node
    #Presentes apenas na primeira + presentes apenas na segunda

    dict_heatmap = {}

    for i in dict_result_edges_1:
        #VERIFICA SE EXISTE NODE IGUAL
        if not dict_equal_nodes[i].empty:
            aux_1 = dict_result_edges_1[i]
            aux_2 = dict_result_edges_2[i]
            
            #IDENTIFICANDO A MAIOR POSIÇÃO DA PROTEÍNA
            maior_pos_protein_1 = dict_result_edges_1[i]['aux_position'].max()
            maior_pos_protein_2 = dict_result_edges_2[i]['aux_position'].max()
            maior = max(maior_pos_protein_1, maior_pos_protein_2)

            #CONTAR O NUMERO DE INTERAÇÕES DE CADA NÓ
            aux_1_series = aux_1['NodeId1'].value_counts()
            aux_1_dataframe = pd.DataFrame({'NodeId1': aux_1_series.index, 'Quant': aux_1_series.values})

            aux_2_series = aux_2['NodeId1'].value_counts()
            aux_2_dataframe = pd.DataFrame({'NodeId1': aux_2_series.index, 'Quant': aux_2_series.values})

            #POPULAR A COLUNA AUX_NODEID E AUX_POSITION
            aux_1_dataframe['aux_nodeId1'] = aux_1_dataframe['NodeId1'].apply(populate_aux_target)
            aux_1_dataframe['aux_position'] = aux_1_dataframe['NodeId1'].apply(populate_aux_position)
            aux_2_dataframe['aux_nodeId1'] = aux_2_dataframe['NodeId1'].apply(populate_aux_target)
            aux_2_dataframe['aux_position'] = aux_2_dataframe['NodeId1'].apply(populate_aux_position)


            
            sum_merge = aux_1_dataframe.merge(aux_2_dataframe, on="aux_nodeId1", how='outer', suffixes=['', '_'], indicator=True)
            
            sum_merge['Comentary'] = ''
            sum_merge['Differences'] = 0
            sum_merge['Color'] = -1

            for w, row in sum_merge.iterrows():
                if row._merge == 'both':
                    sum_merge.loc[w,'Differences'] = int(row.Quant) + int(row.Quant_)
                    sum_merge.loc[w,'Comentary'] = 'Have changes'
                    sum_merge.loc[w,'Color'] = sum_merge.loc[w,'Differences'] + 60 #variação de cor de acordo com a difference
                if row._merge == 'left_only':
                    sum_merge.loc[w,'Differences'] = row.Quant
                    sum_merge.loc[w,'Comentary'] = 'Have changes'
                    #sum_merge.loc[w,'Comentary'] = 'AA not assigned in Protein 2'
                    #sum_merge.loc[w,'Color'] = 30 #apenas presente na proteina 1
                    sum_merge.loc[w,'Color'] = sum_merge.loc[w,'Differences'] + 60 #variação de cor de acordo com a difference
                if row._merge == 'right_only':
                    sum_merge.loc[w,'Differences'] = row.Quant_
                    #sum_merge.loc[w,'Comentary'] = 'AA not assigned in Protein 1'
                    sum_merge.loc[w,'Comentary'] = 'Have changes'
                    sum_merge.loc[w,'aux_position'] = row.aux_position_
                    sum_merge.loc[w,'NodeId1'] = row.NodeId1_
                    #sum_merge.loc[w,'Color'] = 1 #apenas presente na proteina 2
                    sum_merge.loc[w,'Color'] = sum_merge.loc[w,'Differences'] + 60 #variação de cor de acordo com a difference

                            
            
            sum_merge.drop(columns=['Quant', 'NodeId1_', 'Quant_', 'aux_position_', '_merge'], inplace = True) 
            

            #PEGANDO POSIÇÕES DE NÓS SEMELHANTES
            equal_nodes_heatmap = dict_equal_nodes[i].drop(dict_equal_nodes[i].columns[[3,4,5,6,7,8,9,10,11,12,13,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32]], axis=1)
            equal_nodes_heatmap.columns = ['NodeId1', 'Chain', 'Position','aux_nodeId1']
            
            #MERGE COM O DATAFRAME DE NODES IGUAIS PARA PEGAR A POSIÇÃO E NOME DE NÓS PRESENTE NAS DUAS PROTEINAS 
            #E QUE NAO TIVERAM DIFERENÇAS NOS EDGES
            sum_merge = sum_merge.merge(equal_nodes_heatmap, on="aux_nodeId1", how='outer', suffixes=['', '_'], indicator=True)

            for w, row in sum_merge.iterrows():
                if row._merge == 'right_only':
                    sum_merge.loc[w,'NodeId1'] = row.NodeId1_
                    sum_merge.loc[w,'Differences'] = 0
                    sum_merge.loc[w,'Comentary'] = 'No changed'
                    sum_merge.loc[w,'Color'] = 60 #sem diferenças
                    sum_merge.loc[w,'aux_position'] = row.Position

            sum_merge.drop(columns=['NodeId1_', 'Chain', 'Position', '_merge'], inplace = True) 

            #POPULAR ESPAÇOS NAO ASSIMILADOS PELO RING
            cont_df = pd.DataFrame()
            cont_df['aux_position'] = ''
            cont = 1
            
            for w in range(0, int(maior)):
                cont_df.loc[w, 'aux_position'] = cont
                cont+=1

            sum_merge = sum_merge.merge(cont_df, on="aux_position", how='outer', suffixes=['', '_'], indicator=True)

            for w, row in sum_merge.iterrows():
                if row._merge == 'right_only':
                    sum_merge.loc[w,'NodeId1'] = '-'
                    sum_merge.loc[w,'aux_nodeId1'] = '-'
                    sum_merge.loc[w,'Comentary'] = 'AA not assigned in RIN'
                    sum_merge.loc[w,'Differences'] = 0
                    sum_merge.loc[w,'Color'] = 50 #nao assimilado pelo ring

            sum_merge.drop(columns=['_merge'], inplace = True) 

            dict_heatmap[i] = sum_merge.sort_values('aux_position')
            
            
            #VERIFICAÇÃO DE AMINOACIDO EXCLUSIVO DA PROTEINA
            exclusive_1 = {}
            exclusive_2 = {}
            dict_aux_exclusive_1 = dict_result_nodes_1[i]
            dict_aux_exclusive_2 = dict_result_nodes_2[i]
            
            for w, row in dict_aux_exclusive_1.iterrows():
                key_dict = row.Position
                exclusive_1[key_dict] = row.Residue  
                
            for w in exclusive_1:
                index_row = dict_heatmap[i][dict_heatmap[i].aux_position == w].index.tolist()
                dict_heatmap[i].loc[index_row,'Comentary'] = 'AA not assigned in Protein 2'
                dict_heatmap[i].loc[index_row,'Color'] = 30
            
            for w, row in dict_aux_exclusive_2.iterrows():
                key_dict = row.Position
                exclusive_2[key_dict] = row.Residue  
                
            for w in exclusive_2:
                index_row = dict_heatmap[i][dict_heatmap[i].aux_position == w].index.tolist()
                dict_heatmap[i].loc[index_row,'Comentary'] = 'AA not assigned in Protein 1'
                dict_heatmap[i].loc[index_row,'Color'] = 1
                
            #VERIFICAÇÃO DE MUDANÇA DE AMINOACIDO NA POSIÇÃO
            dict_aux_change = dict_result_nodes_change[i]
            dict_aux_change['aux_nodeId'] = dict_aux_change['NodeId'].apply(populate_aux_target)
            
            a = {}
            for w, row in dict_aux_change.iterrows():
                if row.Position not in a:
                    key_dict = row.Position
                    a[key_dict] = row.Residue  
                else:
                    a[key_dict] = a[key_dict] + '-' + str(row.Residue)
            for w in a:
                index_row = dict_heatmap[i][dict_heatmap[i].aux_position == w].index.tolist()
                dict_heatmap[i].loc[index_row,'Comentary'] = 'Amino Changed (' + a[w] + ')'
                dict_heatmap[i].loc[index_row,'Color'] = 800000000
                
                
                
            dict_heatmap[i].drop(columns=['aux_nodeId1', 'aux_position'], inplace = True) 
            
            
            #SAIDA DE ARQUIVO HEATMAP
            file = open(output_directory+"/"+i+'/heatmap-'+i+'.csv', 'a', newline='')
            file.write("group,variable,color,idnode,diff,comentary\n") 

            col = 0
            lin = 1
            #dim_max = 9
            dim_max = int(sqrt(maior))
            pos = 1
            cont = 0

            for w, row in dict_heatmap[i].iterrows():          
                cont+=1
                file.write(str(col) + ',' + str(lin) + ',' + str(row['Color']) + ',' + str(row['NodeId1']) + ',' + str(int(row['Differences'])) + ',' + str(row['Comentary']) + '\n')

                if col <= dim_max:
                    col+=1
                else:
                    lin+=1
                    col = 0
                
                if col == dim_max:
                    lin = cont + 1
                    col = 0
            file.close()


    #ANALISE

    #GRÁFICO 1 - SAÍDA: DATAFRAME ONDE CADA TUPLA É O DEGREE DA POSIÇÃO DE TODAS AS COMPARAÇÕES

    amino_name = ['ALA','ARG','ASN','ASP','CYS','GLU','GLN','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']

    #CRIAÇÃO DE DATAFRAME COM O NOME DE CADA PDB SENDO UMA COLUNA
    df_degree = pd.DataFrame()
    df_degree['Position'] = ''
    for i in dict_protein_nodes_chain.keys():
        df_degree[i]= ''

    #PREENCHIMENTO DA COLUNA POSIÇÃO COM OS NÚMEROS DAS POSIÇÕES SEM REPETIÇOES
    for i in dict_protein_nodes_chain:
        for w, row in dict_protein_nodes_chain[i].iterrows():
            if row["Position"] not in df_degree["Position"].values:
                if row['Residue'] in amino_name:
                    df_degree = df_degree.append({'Position': row["Position"]}, ignore_index = True)

    #ORDENAR COLUNA POSITION E ORDENAR
    df_degree = df_degree.sort_values('Position')
    df_degree.set_index('Position', inplace=True)

    #PREENCHER OS DEGREES DE CADA PDB
    for i in dict_protein_nodes_chain:
        for w, row in dict_protein_nodes_chain[i].iterrows(): 
            if row['Residue'] in amino_name:
                df_degree.loc[row['Position'], i] = row['Degree']

    #REMOVENDO POSIÇÕES BASEADO NA QUANTIDADES DE VALORES NAN
    sum_nan = df_degree.isnull().sum(axis=1)
    df_sum_nan = sum_nan.to_frame()
    for w, row in df_sum_nan.iterrows():
        if (row[0] >= len(df_degree.columns)*80/100) or (abs(row[0] - len(df_degree.columns)) == 1):
            df_degree.drop(w, inplace=True)

    #CRIANDO COLUNAS MEDIA E MEDIANA
    #df_degree['mean'] = df_degree.mean(axis=1)
    #df_degree['median'] = df_degree.median(axis=1)
    df_degree['std'] = df_degree.std(axis = 1, skipna = True)
    df_degree['Position'] = df_degree.index

    file = open(output_directory+ '/df_degree_analysis.txt', 'a', newline='')
    df_degree.to_csv(file, sep="\t", index=False)
    file.close()


    # Data for plotting
    x = df_degree.index
    y = df_degree['std']
    y1 = - df_degree['std']

    # Create two subplots sharing y axis
    fig,ax = plt.subplots()

    #ticks dos pontos dos graficos
    ax.set_xticks(df_degree.index)

    #tamanho do grafico gerado em px
    fig.set_size_inches(130, 10)

    #plot dos dados do dataframe
    ax.plot(x, y, -y, color='xkcd:azure')

    #plot das linhas horizontais
    line1 = plt.axhline(y=0.5, color='xkcd:sand', linestyle='--')
    line2 = plt.axhline(y=-0.5, color='xkcd:sand', linestyle='--')
    line3 = plt.axhline(y=1.5, color='xkcd:peach', linestyle='--')
    line4 = plt.axhline(y=-1.5, color='xkcd:peach', linestyle='--')

    #preenchimento da cor
    ax.fill_between(x, y, y1, color='xkcd:powder blue')


    #plot dos pontos coloridos
    for i in x:
        aux_y = y[i]
        if aux_y >= 1.5:
            plt.plot(i, aux_y, 'ro', color = 'red')
            ax.annotate('  '+str(int(i)), xy=(i, aux_y), xytext=(1, 0), textcoords="offset points", ha='left', va='center', color='red', clip_on=True)
        elif aux_y >= 0.5 and aux_y < 1.5:
            plt.plot(i, aux_y, 'ro', color = 'xkcd:orangey brown')
            ax.annotate('  '+str(int(i)), xy=(i, aux_y), xytext=(1, 0), textcoords="offset points", ha='left', va='center', color='xkcd:orangey brown', clip_on=True)
        else:
            plt.plot(i, aux_y, 'ro', color = 'xkcd:cerulean')
            ax.annotate('  '+str(int(i)), xy=(i, aux_y), xytext=(1, 0), textcoords="offset points", ha='left', va='center', color='xkcd:cerulean', clip_on=True)
    for i in x:
        aux_y = y1[i]
        if aux_y <= -1.5:
            plt.plot(i, aux_y, 'ro', color = 'red')
            ax.annotate('  '+str(int(i)), xy=(i, aux_y), xytext=(1, 0), textcoords="offset points", ha='left', va='center', color='red', clip_on=True)
        elif aux_y <= -0.5 and aux_y > -1.5:
            plt.plot(i, aux_y, 'ro', color = 'xkcd:orangey brown')
            ax.annotate('  '+str(int(i)), xy=(i, aux_y), xytext=(1, 0), textcoords="offset points", ha='left', va='center', color='xkcd:orangey brown', clip_on=True)
        else:
            plt.plot(i, aux_y, 'ro', color = 'xkcd:cerulean')
            ax.annotate('  '+str(int(i)), xy=(i, aux_y), xytext=(1, 0), textcoords="offset points", ha='left', va='center', color='xkcd:cerulean', clip_on=True)
            
    #titulo
    ax.set(title='Titulo', ylabel='Standart Deviation')

    #ajustas as bordas do gráfico para tocar nos eixos
    plt.axis(xmin=df_degree.index.min() - 0.2, xmax=df_degree.index.max())

    #adiciona o grid
    ax.grid()

    plt.savefig(output_directory+'/plot.svg')


    directory = "../" + directory
    r_script_path = os.path.join(project_path, "scripts_corins", "code.R")
    pipe = subprocess.Popen(["Rscript", r_script_path, directory])
    pipe.wait()
    # #subprocess.call ("Rscript --vanilla ./code.R", shell=True)



if (__name__ == "__main__"):
    # nem aqui tbm n ta sendo
    print("äaaaaaa   ", sys.argv[1])
    script(sys.argv[1], sys.argv[2])


#script("../../upload_ring_files/")

