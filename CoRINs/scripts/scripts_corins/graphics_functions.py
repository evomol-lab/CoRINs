from math import sqrt

#retonar dicionarios e cadeias do dicionarios
#nome do arquivo
def criarDicionarioNode(nomeArq):

    file = open(nomeArq, "r")

    file.readline() # ignora a primeira linha do arquivo
    t = file.readlines()

    dict = {}
    chains = []
    
    #
    for linha in t:
        linha = linha.split("\t")
 
        if (linha[1] not in dict):
            dict[linha[1]] = [[int(linha[2]), linha[0], int(linha[5])]]
            chains.append(linha[1])            
        else:
            dict[linha[1]].append([int(linha[2]), linha[0], int(linha[5])])

    for i in dict:
        dict[i].sort()
	
    return dict, chains

def criarDicionarioEdge(nomeArq, param1, param2):

    file = open(nomeArq, "r")    
    
    dict = {}    
    
    # dict[head_file] = {} <- um Dicionarios de dicionario 

    for linha in file.readlines():
        linha = linha.split('\t')
        
        if ("NodeId1" in linha[0]):
            
            a = linha[0].split(' ')[-1]
            b = linha[-1].split(' ')[1]

            head_file = param1 + "_" + a + "-" + b + "_" + param2
            dict[head_file] = {}    
        else:
            if (linha[0] in dict[head_file]):
                dict[head_file][linha[0]] += 1
            else :
                dict[head_file][linha[0]] = 1

            if (linha[2] in dict[head_file]):
                dict[head_file][linha[2]] += 1
            else :
                dict[head_file][linha[2]] = 1 

    # se necessario tem que ordenar aqui 
    return dict
    
# cria grafico degree_heatmap
def degree_heatmap (namefile, dict, path_result):
    
    
    for chain in dict:
        
        name_path = path_result + "/heat_degree_"+namefile+"_"+chain+".csv"
        file2 = open(name_path, "w")
        file2.write("group,variable,idnode,degree,color\n")       
        dict[chain] = sorted(dict[chain], reverse=True)
 
        linha = 0
        col = 1
           
        

        for i in range(1, int(dict[chain][0][0]) + 1):
            nome = "-"
            degree = 0 
            color = 75           

            for j in dict[chain]:
                if i > j[0]:
                    break
                elif i == j[0]:
                    nome = j[1]
                    degree = j[2]
                    color = int(degree * 1) + 75

            file2.write(str(col) + "," + str(linha) +  "," + nome + "," + str(degree) + ","+ str(color) +"\n")

            col += 1

            if col == int(sqrt(int(dict[chain][0][0]))):
                col = 1
                linha += 1
           
#retonar a maior sequencia das cadeias
def compare_chain(dict1, dict2):
	
	x = dict1[0][0]
	y = dict2[0][0]
	
	return max(x, y)

# cria grafico diff_degree
def count_diff_degree (diretorio, path_result):
    file = open(diretorio, "r")

    array = []

    for line in file.readlines():
        line = line.split() 

        if len(line) == 7:
            nome = path_result + "/degree-"+ line[0] +"x" + line[-1] + "-" + line[2] + line[4] +".tsv"
            file = open( nome , "w")  
            file.write("nodes\tdifferences\n")  
            x = ""   
            y = ""
            array = []

        elif (len(line) == 6 and line[0] != "NodeId"):

                if x == "":
                    x = line
                else :
                    y = line                              
                    #file.write(x[0] + " " +str(abs(int(x[3]) - int(y[3]))) + "\n")
                    array.append([abs(int(x[3]) - int(y[3])), x[0]])
                    x = ""
                    y = ""
        elif len(line) == 3:
            if len(array) != 0:
                array = sorted(array, reverse=True)                
                for i in range(0, 10):
                    file.write(array[i][1] + "\t" +str(array[i][0]) + "\n")

# cria grafico diff_interection
def count_diff_interaction (diretorio, nome, path_result) :

    file = open(diretorio, "r")

    dict = {}
    
    keys = []

    for line in file.readlines():
        line = line.split()

        if line[0] == "NodeId1":
            chave = line[2]+line[4]            
            keys.append(chave)
            dict[chave] = []
        else :
            dict[chave].append([line[1]])

    for k in keys:        
        dict2 = {}
        
        for line in dict[k]:            
            for i in line:
                if i not in dict2:
                    dict2[i] = 1
                else :
                    dict2[i] += 1
        
        path = path_result + "/"+nome+"-"+k+".tsv"
        file2 = open(path, "w")
        file2.write("interactions\tlosses\n")

        array = []
        
        for i in dict2:
            array.append([dict2[i], i])

        array = sorted(array, reverse=True)

        for i in array:
            file2.write(i[1] + "\t" + str(i[0]) + "\n")
                
# No exemplo \Entrada\Eden uma amino tem a mesma posição e nome diferentes 
# o caso acontece no 1aki A x A 1at6 (A:101:_:ASP) <3
def heatmap (param, array_dict_arq, array_key_arq, dict_arq_edge, directory):
    
    #cont = 0 (apenas para teste)

 

    for i in range(0, len(param)):
        
        for j in range(0, len(param)):

            if i != j :
                
                #print("\n" + str(cont) + "\n")  apenas para teste
                for l in array_key_arq[i]:
                    for k in array_key_arq[j]:
                        #cont += 1 (apenas para teste)
                        #print(cont) (apenas para teste)
                        #print(param[i] + " " + l+" x "+ k + " " + param[j]) (apenas para teste)                        
                        
                        # encontra o maior valor entre as cadeias
                        maior_da_cadeia = compare_chain(array_dict_arq[i][l], array_dict_arq[j][k])
                        
                        #
                        path = directory + 'heatmap-'+param[i]+'x'+param[j]+'-'+l+k+'.csv'
                        file = open(path, 'w')

                        file.write('group,variable,color,idnode,diff,comentary\n')

                        array_diff_node_1 = {}
                        array_diff_node_2 = {}
                        
                        # Verificar a ausencia e presença de nós
                        # Presente na primeira e ausente na segunda
                        for x in array_dict_arq[i][l]:
                            a = True
                            aux1 = x[1][1::]
                            for y in array_dict_arq[j][k]:
                                aux2 = y[1][1::]
                                

                                if aux1 == aux2:
                                    a = False
                                    break
                            if a:
                                array_diff_node_2[x[0]] = x[1]

                        

                        # Presente na segunda e ausente na primeira
                        for x in array_dict_arq[j][k]:
                            a = True
                            aux1 = x[1][1::]
                            for y in array_dict_arq[i][l]:
                                aux2 = y[1][1::]
                                if aux1 == aux2:
                                    a = False
                                    break
                            if a:
                                array_diff_node_1[x[0]] = x[1]

                        
                        
                                                       
                                                           
                        # contagem de diferença correspondente aos valores das cores	
                        
                        coluna = 1
                        linha = 0
                        head_file = param[i] + "_" + l + "-" + k + "_" + param[j]

                        #print(dict_arq_edge[param[i]+"-"+param[j]][head_file]) (apenas para teste)


                        for cont in range(1, maior_da_cadeia + 1):
                            nome = '-'

                            # procurando o nome do node no dicionario 1
                            for a in array_dict_arq[i][l]:
                                if cont == a[0]:
                                    nome = a[1]

                            # procurando o nome do node no dicionario 2
                            if nome == '-':
                                for a in array_dict_arq[j][k]:
                                    if cont == a[0]:
                                        nome = a[1]                        
                            

                            cor = 60
                            diff = 0

                            if (nome != '-') and (nome in dict_arq_edge[param[i]+"-"+param[j]][head_file]):                                
                                diff = dict_arq_edge[param[i]+"-"+param[j]][head_file][nome]
                                cor += (dict_arq_edge[param[i]+"-"+param[j]][head_file][nome] * 2.5)

                            if cont in array_diff_node_1:   #  Não tem na primeira e tem na segunda
                                file.write(str(coluna) + "," + str(linha) + "," + str(1) +",\"" + nome +"\"," + str(diff) + "\"," + "AA not assigned in Protein 1" +"\n")
                            elif cont in array_diff_node_2: # Tem na primeira e nao tem na segunda
                                file.write(str(coluna) + "," + str(linha) + "," + str(30) +",\"" + nome +"\"," + str(diff) +"\"," + "AA not assigned in Protein 2" +"\n")
                            elif nome == "-":               # nao tem na primeira e nao tem na segunda
                                file.write(str(coluna) + "," + str(linha) + "," + str(50) +",\"" + nome +"\"," + "-" + "\"," + "AA not assigned in RIN" +"\n")
                            else:                           # tem na primeira e tem na segunda
                                file.write(str(coluna) + "," + str(linha) + "," + str(cor) +",\"" + nome +"\"," + str(diff) + "\"," + " " +"\n")

                            coluna += 1
                            if coluna == int(sqrt(maior_da_cadeia)):
                                linha = cont + 1
                                coluna = 1

                            

                            





            