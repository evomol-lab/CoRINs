# CoRINs (Comparator of Residue Interaction Networks)
O CoRINs é uma ferramenta Web de comparação de múltiplos arquivos provenientes de redes de interação de resíduos do software RING 2.0 (https://protein.bio.unipd.it/ring/) gerando resultados de fácil compreensão através de tabelas, gráficos interativos e arquivos de texto. Desenvolvida em linguagem Python com bibliotecas de alta performance como Pandas (https://pandas.pydata.org/), Numpy (https://numpy.org) e Matplotlib (https://matplotlib.org/), juntamente com os frameworks de desenvolvimento Web Django (www.djangoproject.com/), Bootstrap (https://getbootstrap.com/), D3.js (https://d3js.org/) e um módulo de geração de parâmetros de grafos (betweenness e coeficiente de clusterização) desenvolvido em linguagem R (https://www.r-project.org/).

### Autores: 
 - Felipe Vieira da Fonseca (Author, Data Scientist)
 - Thiago Dantas Soares (Web Developer, Technical Support)
 - Diego Morais (R Developer <betweenness e coeficiente de clusterização script>)
 - João Paulo Matos Santos Lima (Coordinator, Author, Maintainer).

### Citation: 
https://www.biorxiv.org/------------

# How does it work?
A entrada dos dados do CoRINs consiste nos arquivos nodes.txt e edge.txt ou sua forma compactada gerada pela versão web do RING 2.0 (https://protein.bio.unipd.it/ring/), de cada RIN a ser comparada. Os nodes e edges são convertidos em dataframes, separados por cadeias proteicas e comparados um a um de acordo com a posição dos aminoácidos e das diferentes cadeias das proteínas. Assim são obtidos os resultados de acordo com a comparação dois a dois de cada cadeia para cada proteína, incluindo a comparação entre as cadeias de uma mesma proteína.
A partir da tela inicial do CoRINs, é possível efetuar o carregamento de todas as RINs a serem comparadas usando o botão Load, que abrirá um gerenciador de arquivos de acordo com o sistema operacional utilizado. Após todas as proteínas serem carregadas no CoRINs, os resultados são gerados e organizados por comparação e os dados são mostrados em gráficos e tabelas, conforme o tipo de análise. 

# Interpreting Results
O CoRINs é capaz de mostrar algumas informações úteis para a análise das redes, de forma global, ou seja, uma comparação entre todos as RINs carregadas e local, comparando as RINs cadeia a cadeia.

* Standard deviation of degree of all comparisons

  ![alt text](https://github.com/ffvieira18/teste/blob/master/Images/grafico1-standard-deviation-of-degree-of-all-comparisons.png)

  A análise global consiste em um gráfico que representa a variação na formação de interações de todos os resíduos. O gráfico é construído a partir do desvio padrão de todos os valores de degrees dos aminoácidos equivalentes em posição, respeitando a condição da existência de 80% e ao menos dois degrees na posição referida.
Este gráfico apresenta como eixo x, a posição do aminoácido e como eixo y, o desvio padrão com o seu resultado refletido negativamente para melhor visualização da variação do degree na posição referida. Também foram traçadas duas linhas guias, nos valores de desvio de padrão de 0.5 e 1.5, onde foram estabelecidos os critérios das cores dos pontos (Tabela). 

  | Intervalo        | Cor          | 
  | ---------------- |:------------:|
  | >= 0.0 e <= 0.5  | Azul         |
  | > 0.5 e < 1.5    | Marrom       | 
  | >= 1.5           | Vermelho     |   

  Desvios padrões em torno do valor de degree médio com variação acima de 0.5 já indicam resíduos que variam consideravelmente o seu número de conexões. Valores de desvio padrão maiores que 1.5 indicam prováveis posições de variação conformacional.

* Total differences in interactions

  ![alt text](https://github.com/ffvieira18/teste/blob/master/Images/grafico2-Total-differences-in-interactions.png)

  O gráfico do total de diferenças das interações, consiste no número total de interações químicas que foram perdidas entre a comparação de duas RINs, portanto uma análise local. São mostradas apenas as interações estabelecidas pelo software RING 2.0, sendo as mais comuns: VDW, HBOND, PIPISTACK, PICATION e IONIC. 
O eixo x apresenta o tipo de interação e o eixo y, o número total de interações perdidas. Este gráfico é gerado a partir da biblioteca D3.js e possui interatividade com o cursor do mouse, mudando a cor do box de azul para vermelho.

* High number degree

  ![alt text](https://github.com/ffvieira18/teste/blob/master/Images/grafico3-High-number-degree.png)
  
  Para melhor visualização dos valores de degree, o CoRINs gera dois gráficos no formato de heatmap, um para cada RIN da comparação, com o objetivo de trazer uma visualização ampla dos degrees de cada aminoácido na RIN e orientar o usuário sobre quais os resíduos de degrees mais altos ou com diferenças mais altas em seus valores.
  Ao passar com o mouse no heatmap, é mostrado o Id do aminoácido, composto pela cadeia, posição e a sigla do nome do aminoácido, como também o seu valor de degree. A escala de coloração desse gráfico depende do valor associado ao degree daquela posição (Imagem), do mais claro, próximo ao branco, a inexistência de aminoácido na posição, até uma coloração mais escura, até o vermelho, para os degrees mais altos.

  <img src="https://raw.githubusercontent.com/d3/d3-scale-chromatic/master/img/OrRd.png" width="100%" height="40" alt="OrRd">
  
  Com o objetivo de facilitar a visualização, o eixo x é representado pelo número da primeira posição do aminoácido da linha e o eixo y por uma sequência de números de 0 até um determinado número, equivalente a raiz quadrada do total de aminoácidos. Como exemplo, para localizar o aminoácido da posição 124, podemos nos orientar pela linha de um número de posição próximo, como a linha de valor 121 e realizar a soma com um número equivalente do eixo y em que o resultado seja 124, no caso o valor 3.

* High number degree difference

  ![alt text](https://github.com/ffvieira18/teste/blob/master/Images/grafico4-high-numbers-of-degree-differences.png)

  Pode ser interessante para o usuário observar rapidamente quais os aminoácidos que mais variaram em degree. Para isso, o CoRINs gera um gráfico utilizando a biblioteca D3.js com as 10 maiores diferenças nos degrees da rede (Figura 12). Sendo o eixo x o Id do aminoácido e o eixo y o valor da diferença. Por exemplo, o aminoácido leucina da posição 168 da cadeia A, possui degree de valor 2 na RIN da proteína de PDB 1B55 e valor 6 na proteína de PDB 2z0P, portanto, possui uma diferença de 4 no degree deste aminoácido.   

* Interactive Tables (Differences of Degrees, Differences in Aminoacid, Differences of Nodes, Differences of Edges, Graph parameters)

  ![alt text](https://github.com/ffvieira18/teste/blob/master/Images/grafico5-differences-of-degrees.png)

  Para melhor representação dos resultados, foram implementadas tabelas interativas com a biblioteca datatables (https://datatables.net/). Assim é possível verificar facilmente quais aminoácidos variam em Degree, troca de aminácidos, ausência ou presença de aminoácidos e coeficientes de clusterização e betweenness de cada aminoácido da rede. As tabelas também apresentam um campo Search para pesquisar por qualquer caractere ou palavra que esteja presente no Id do aminoácido e como resultado será mostrada apenas a lista de aminoácidos contendo a palavra pesquisada; O campo Show entries, onde é possível alterar o número de linhas do resultado e os botões para alterar a páginas do resultado.

* Visualization of Difference RIN

  ![alt text](https://github.com/ffvieira18/teste/blob/master/Images/grafico6-visualization-of-differences-RIN.png)

  Para uma visualização geral da comparação, foi gerado um gráfico semelhante ao gráfico da visualização geral dos valores de degrees. Neste gráfico é possível observar por diferentes cores, as posições que se encontram com aminoácidos ausentes na cadeia, a presença exclusiva de um aminoácido em apenas uma das RINs, troca de aminoácidos em uma determinada posição e quantidade de diferenças nas interações químicas (Tabela). Devido ao número de posições, ou pelo fato dos ligantes e heteroátomos apresentarem posições muito elevadas. Assim, é sempre aconselhável marcar a opção para retirar os ligantes e heteroátomos pelo software RING 2.0. 

  | Descrição       | Cor          | 
  | ---------------- |:------------:|
  | Aminoácido ausente na posição nas duas RINs  | Cinza         |
  | Quantidade de diferenças nas interações químicas do aminoácido daquela posição    | Escala de laranja-rosado       | 
  | Aminoácido presente apenas na RIN 1           | Azul Escuro     |   
  | Aminoácido presente apenas na RIN 2           | Azul Claro   |  
  | Posição em que ocorreu troca de aminoácido           | Vermelho     |  

# Installing on Windows

Install Docker
```
<https://www.docker.com/>
```
Open CMD, enter the project folder and execute docker build command
```
<docker-compose build>
```
Type the command to run the application
```
<docker-compose up>
```
# Installing on Linux

Install Docker
```
<sudo apt-get install -y docker.io>
```
Install Docker-Compose
```
<sudo curl -L "https://github.com/docker/compose/releases/download/1.26.0/docker-compose-$(uname -s)-$(uname -m)" -o /usr/local/bin/docker-compose>
<sudo chmod +x /usr/local/bin/docker-compose>
<sudo usermod -aG docker "$USER">
```
Open CMD, enter the project folder and execute docker build command
```
<docker-compose build>
```
Type the command to run the application
```
<docker-compose up>
```
