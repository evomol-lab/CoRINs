""" from graphics_functions import *

dict_arq1, chaves_arq1 = criarDicionario("../../upload_ring_files/1b55_nodes.txt")
dict_arq2, chaves_arq2 = criarDicionario("../../upload_ring_files/2z0p_nodes.txt")

#
degree_heatmap("../../upload_ring_files/1b55_nodes.txt", dict_arq1)

# os maiores das cadeias
cadeia = comparaCadeia(dict_arq1, dict_arq2)

# Cria dicionarios dos nodes
file1 = "../../result/result_corins/dif_nodes_2z0p_1b55.txt-"
dict_node_diff1 = diferencaAmino(chaves_arq1, chaves_arq2, file1)

file2 = "../../result/result_corins/dif_nodes_1b55_2z0p.txt-"
dict_node_diff2 = diferencaAmino(chaves_arq2, chaves_arq1, file2)

# Cria dicionarios dos Edges
file1 = "../../result/result_corins/dif_edges_1b55_2z0p.txt"
dict_edge_diff1 = count_diff_edge(file1)

file2 = "../../result/result_corins/dif_edges_2z0p_1b55.txt"
dict_edge_diff2 = count_diff_edge(file2)

# Cria os arquivo do heatmap
nome1 = "1b55x2z0p"
output_svg(chaves_arq1, chaves_arq2, dict_arq1, dict_arq2, cadeia, dict_node_diff1, dict_node_diff2, dict_edge_diff1,nome1)

nome2 = "2z0px1b55"
output_svg(chaves_arq2, chaves_arq1, dict_arq2, dict_arq1, cadeia, dict_node_diff2, dict_node_diff1, dict_edge_diff2,nome2)

# criar um dicionarios as diferenças dos degree de cada node
file1 = "../../result/result_corins/dif_degree_1b55_2z0p.txt"
count_diff_degree(file1)

file2 = "../../result/result_corins/dif_degree_2z0p_1b55.txt"
count_diff_degree(file2)


# criar um dicionarios as diferenças das interaction de cada edge
file1 = "../../result/result_corins/dif_edges_1b55_2z0p.txt"
count_diff_interaction(file1, nome1)

 """