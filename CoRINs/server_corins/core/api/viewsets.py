from rest_framework.views import APIView
from rest_framework.viewsets import ModelViewSet
from rest_framework.decorators import action
from rest_framework.response import Response
from core.models import Arquivo
from .serializers import CoreSerializers
import os, shutil, re, json
from glob import glob
import subprocess
import pandas as pd
from subprocess import PIPE, Popen
from django.conf import settings
import io

class CoreViewSet(ModelViewSet):
    queryset = Arquivo.objects.all()
    serializer_class = CoreSerializers

    def create(self, request, *args, **kwargs):

        
        if not request.session.exists(request.session.session_key):
            request.session.create()

        pasta_usuario = request.session.session_key
        pasta_usario_projeto = '../upload_ring_files/'+pasta_usuario+'/'
        arqvs = []
        for arqv in request.FILES.getlist('files'):
            instance = Arquivo(files = arqv)
            instance.save()
            arqvs.append(arqv.name)
        
        os.mkdir(pasta_usario_projeto)

        for arq in arqvs:
            
            shutil.move('./up_dj_ring_files/'+arq,  pasta_usario_projeto)
        
        print("bbbbbb", pasta_usario_projeto)
        felipe_script_path = os.path.join(settings.BASE_DIR, "scripts_corins", "Felipe.py")
        pipe = subprocess.Popen(["python", felipe_script_path, pasta_usario_projeto, settings.BASE_DIR])
        pipe.wait()

    

        #report_dif_nodes_
        regex = r"^report_dif_nodes_.*"
        x=[]
        # for arquivos in os.walk(pasta_usario_projeto):
            # matches = re.findall(regex, arquivos, re.MULTILINE)
            # if len(matches) > 0:
            #     x.append(matches[0])
        listagem = []
        arquivos =subfolders = [f.path for f in os.scandir(pasta_usario_projeto) if f.is_dir() ] 
        for arq_str in arquivos:
            arq = arq_str.split('/')[-1]
            listagem.append(arq)
        print(listagem)
        # print(x)
        # for arq in x:

        #     listagem.append(arq.split("_")[-1][:-4])
        
        # print(listagem)
         
        
        response = {
          'nome_arquivos' : arqvs,
          'pasta': pasta_usuario,
          'arquivos': listagem
        }
        
        print(response)
        
        return Response( response )

    @action(methods=['GET'], detail=False)
    def pegarchains(self, request, pk=None):
        
        pasta = request.GET['pasta']
        nome_prot = request.GET['prot']
        
        print(nome_prot)

        pasta_usario_projeto = '../upload_ring_files/'+pasta+'/'+nome_prot+'_nodes.txt'
        dict_chain = {}
        chain = []
        file = open(pasta_usario_projeto,'r')
        file.readline()
        for linha in file.readlines():
            linha = linha.split("\t")

            if (linha[1] not in dict_chain):
                dict_chain[linha[1]] = 0
                chain.append(linha[1])
        print(chain)
        return Response({'chain':chain})

    @action(methods=['POST'], detail=False)
    def degree(self, request, pk=None):
        
        folder = request.POST['prot1']
        # segunda_prot = request.POST['prot2']

        # chain_1 = request.POST['chain1']
        # chain_2 = request.POST['chain2']
        # print(request.POST)
        pasta = request.POST['pasta']

        # key_to_return = primeira_prot+'-'+chain_1+'x'+chain_2+'-'+segunda_prot

        directory = '../upload_ring_files/'+pasta+'/'+folder+'/report_dif_degree_'+folder+'.txt'
        # print('777777777777777777777777777777777777777777'+directory)
        # file = open(directory, 'r')


        # #1b55 - A x A - 2z0p
        # #NodeId			Chain		Position	Degree		BfactorCA		pdbFileName
        # #A:5:_:ILE		A			5			5			40.500			1b55.pdb#5.A
        # #Total: 102 Aminoacids

        # dict = {}
        # for line in file.readlines():
        #     line = line.split()
        #     if len(line) == 7:
        #         key = line[0] + line[1] + line[2] + line[3] + line[4] + line[5] + line[6]
        #         dict[key] = []

        #     if len(line) != 0 and len(line) < 7:
        #         if line[0] != 'NodeId' and line[0] != 'Total:':
        #             linha = line[0] + '\t' + line[3] + '\t' + line[5]
        #             dict[key].append(linha)

        data = pd.read_csv(directory, sep="\t", header= 0)
        if 'Chain' in data.columns:
            data.drop(columns=['Chain', 'Position', 'Residue', 'Dssp', 'Bfactor_CA', 'x', 'y', 'z', 'Rapdf', 'Tap', 'Accessibility'], inplace = True) 
        # print(data)
        # result = json.dumps(data)

        return Response({'degree':data})
    
    @action(methods=['POST'], detail=False)
    def rcode1(self, request, pk=None):
        
        folder = request.POST['prot1'].split('-')[0][:-3]
        pasta = request.POST['pasta']

        
        directory = '../upload_ring_files/'+pasta+'/ResultadoNos_'+folder+'_edges.txt'
        
        
        data = pd.read_csv(directory, sep="\t", header= 0)
        # data = open('../upload_ring_files/ewnyrcx7kuh4buo88aylr303ly8dgfz2/ResultadoNos_3OG7_edges.txt','r')
        

        return Response({'rcode1':data, 'name': folder})

    @action(methods=['POST'], detail=False)
    def rcode2(self, request, pk=None):
        
        folder = request.POST['prot1'].split('-')[1][:-3]
        pasta = request.POST['pasta']

        
        directory = '../upload_ring_files/'+pasta+'/ResultadoNos_'+folder+'_edges.txt'
        
        
        data = pd.read_csv(directory, sep="\t", header= 0)
        # data = open('../upload_ring_files/ewnyrcx7kuh4buo88aylr303ly8dgfz2/ResultadoNos_3OG7_edges.txt','r')
        

        return Response({'rcode2':data, 'name': folder})
    
    @action(methods=['POST'], detail=False)
    def edgestables(self, request, pk=None):
        
        folder = request.POST['prot1']
        pasta = request.POST['pasta']

        
        directory = '../upload_ring_files/'+pasta+'/'+folder
        file = open(directory+'/report_dif_edges_'+folder+'.txt', 'r') 

        name = file.readline()
        protein_name_1 = name.split('-')[0]
        protein_name_2 = name.split('-')[1]

        #PROTEINA 1
        aux_file = open(directory+'aux_file1.txt', 'w')
        line = file.readline()
        line = file.readline()
        line = file.readline()
        while line:
            line = file.readline()
            if protein_name_2 in line:
                break
            aux_file.write(line)    
        aux_file.close() 
        df_protein_edges_1 = pd.read_csv(directory+'aux_file1.txt', sep="\t", header= 0)
        if 'aux_chain_1' in df_protein_edges_1.columns:
            df_protein_edges_1.drop(['aux_chain_1', 'aux_chain_2', 'Quant', 'aux_position'], inplace=True, axis=1)

        #PROTEINA 2
        aux_file = open(directory+'aux_file2.txt', 'w')
        while line:
            line = file.readline()
            if protein_name_1 in line:
                break
            aux_file.write(line)
        aux_file.close() 
        df_protein_edges_2 = pd.read_csv(directory+'aux_file2.txt', sep="\t", header= 0)
        if 'aux_chain_1' in df_protein_edges_2.columns:
            df_protein_edges_2.drop(['aux_chain_1', 'aux_chain_2', 'Quant', 'aux_position'], inplace=True, axis=1)
        
        
        # data = open('../upload_ring_files/ewnyrcx7kuh4buo88aylr303ly8dgfz2/ResultadoNos_3OG7_edges.txt','r')
        

        return Response({'edge1':df_protein_edges_1,'edge2':df_protein_edges_2, 'name1': protein_name_1, 'name2': protein_name_2[:-1]})
    
    @action(methods=['POST'], detail=False)
    def nodestables(self, request, pk=None):
        
        folder = request.POST['prot1']
        pasta = request.POST['pasta']

        
        directory = '../upload_ring_files/'+pasta+'/'+folder
        file = open(directory+'/report_dif_nodes_'+folder+'.txt', 'r')  

        name = file.readline()
        protein_name_1 = name.split('-')[0]
        protein_name_2 = name.split('-')[1]
        
        #PROTEINA 1
        aux_file = open(directory+'aux_file3.txt', 'w')
        line = file.readline()
        line = file.readline()
        line = file.readline()
        while line:
            line = file.readline()
            if protein_name_2 in line:
                break
            aux_file.write(line)    
        aux_file.close() 
        df_protein_nodes_1 = pd.read_csv(directory+'aux_file3.txt', sep="\t", header= 0)
        if 'Chain' in df_protein_nodes_1.columns:
            df_protein_nodes_1.drop(['Chain', 'Position','Residue', 'Dssp', 'x', 'y', 'z', 'Rapdf', 'Tap', 'Accessibility'], inplace=True, axis=1)
        # print(1)
        # print(df_protein_nodes_1)

        #PROTEINA 2
        aux_file = open(directory+'aux_file4.txt', 'w')
        while line:
            line = file.readline()
            if protein_name_1 in line:
                break
            aux_file.write(line)
        aux_file.close() 
        df_protein_nodes_2 = pd.read_csv(directory+'aux_file4.txt', sep="\t", header= 0)
        if 'Chain' in df_protein_nodes_2.columns:
            df_protein_nodes_2.drop(['Chain', 'Position','Residue', 'Dssp', 'x', 'y', 'z', 'Rapdf', 'Tap', 'Accessibility'], inplace=True, axis=1)
        # print(2)
        # print(df_protein_nodes_2)
        #AMINOCHANGED
        aux_file = open(directory+'aux_file5.txt', 'w')
        while line:
            line = file.readline()
            if line == '':
                break
            aux_file.write(line)
        aux_file.close() 
        df_protein_nodes_3 = pd.read_csv(directory+'aux_file5.txt', sep="\t", header= 0)
        # if 'Chain' in df_protein_nodes_3.columns:
        df_protein_nodes_3.drop(['Chain', 'Position','Residue', 'Dssp', 'x', 'y', 'z', 'Rapdf', 'Tap', 'Accessibility'], inplace=True, axis=1)
        # print(df_protein_nodes_3)
        
        # data = open('../upload_ring_files/ewnyrcx7kuh4buo88aylr303ly8dgfz2/ResultadoNos_3OG7_edges.txt','r')
        

        return Response({'node1':df_protein_nodes_1,'node2':df_protein_nodes_2, 'node3':df_protein_nodes_3 ,'name1': protein_name_1, 'name2':protein_name_2[:-1]})
    
    @action(methods=['POST'], detail=False)
    def nodes(self, request, pk=None):
        #tabela
        primeira_prot = request.POST['prot1']
        segunda_prot = request.POST['prot2']

        chain_1 = request.POST['chain1']
        chain_2 = request.POST['chain2']

        pasta = request.POST['pasta']

        key_to_return = primeira_prot+'-'+chain_1+'x'+chain_2+'-'+segunda_prot

        directory = '../upload_ring_files/'+pasta+'/dif_nodes_'+primeira_prot+'_'+segunda_prot+'.txt'
        file = open(directory, 'r')


        #1b55 - A x A - 2z0p
        #NodeId			Chain		Position	Degree		BfactorCA		pdbFileName
        #A:5:_:ILE		A			5			5			40.500			1b55.pdb#5.A
        #Total: 102 Aminoacids

        dict = {}
        file = open(directory, 'r')

        for line in file.readlines():
            line = line.split()
            if len(line) == 7:
                key = line[0] + line[1] + line[2] + line[3] + line[4] + line[5] + line[6]
                dict[key] = []

            if len(line) != 0 and len(line) < 7:
                if line[0] != 'NodeId' and line[0] != 'Total:' and line[0] != 'Similaridade:':
                    linha = line[0] + '\t' + line[3] + '\t' + line[4]
                    dict[key].append(linha)
        
        return Response({'degree':dict[key_to_return]})

    @action(methods=['POST'], detail=False)
    def edges(self, request, pk=None):
        #tabela
        primeira_prot = request.POST['prot1']
        segunda_prot = request.POST['prot2']

        chain_1 = request.POST['chain1']
        chain_2 = request.POST['chain2']

        pasta = request.POST['pasta']

        key_to_return = chain_1+'x'+chain_2

        directory = '../upload_ring_files/'+pasta+'/dif_edges_'+primeira_prot+'_'+segunda_prot+'.txt'
        file = open(directory, 'r')


        #1b55 - A x A - 2z0p
        #NodeId			Chain		Position	Degree		BfactorCA		pdbFileName
        #A:5:_:ILE		A			5			5			40.500			1b55.pdb#5.A
        #Total: 102 Aminoacids

        dict = {}
        file = open(directory, 'r')

        for line in file.readlines():
            line = line.split()
            if line[0] == "NodeId1":
                key = line[2] + "x" +line[4]
                dict[key] = []

            if line[0] != 'NodeId1':
                linha = line[0] + '\t' + line[1] + '\t' + line[2]
                dict[key].append(linha)
        
        return Response({'degree':dict[key_to_return]})


    @action(methods=['POST'], detail=False)
    def aaa(self, request, pk=None):
        return Response({'chain':1})
