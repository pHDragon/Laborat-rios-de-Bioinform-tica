# Laboratórios de Bioinformática
<br />

## *Staphylococcus aureus* N315
<br />

### Mestrado em Bioinformática

<img src="https://upload.wikimedia.org/wikipedia/commons/9/93/EEUMLOGO.png" width="150" height="136" />
<br />

&nbsp; &nbsp; &nbsp; **Grupo 6:**
* Andreia Rodrigues (PG12573)
* Pedro Moreira (PG38277)
* Pedro Araújo (PG37044)
<br />

## Introdução

*Staphylococcus aureus* é uma bactéria gram-positiva, com forma arredondada, pertencente ao filo Firmicutes e é um membro frequente da microbiota do corpo humano. Normalmente é um organismo comensal, podendo ser um agente patogénico, envolvido em diversas doenças, como infeções da pele, pneumonia, meningite, osteomielite, sépsia, endocardite, entre outros. As estirpes patogénicas promovem infeções através da produção de fatores de virulência. O aparecimento de estirpes de *S. aureus* resistentes a antibióticos é um problema a nível mundial.


O tratamento preferencial para infeções de *S. aureus* é a penicilina, que inibe a síntese das ligações cruzadas entre os aminoácidos do peptidoglicano. No entanto, em muitos países a resistência a penicilina é extremamente comum. A alternativa é frequentemente um antibiótico de β-lactama resistente a penicilase. Os antibióticos β-lactâmicos resistentes, ainda são utilizados como primeira linha de tratamento. A meticilina foi o primeiro antibiótico a ser usado, sendo introduzido em 1959, mas apenas dois anos depois, o primeiro caso de *S. aureus* resistente à meticilina (MRSA).


Devido a esta resistência, foram desenvolvidos antibióticos não β-lactâmicos, como a clindamicina e trimetoprim/sulfametoxazole. Resistência a estes antibióticos também apareceu, que levou ao uso de novos antibióticos abrangentes para bactérias gram-positivas. Antibióticos usados em tratamento de infeções devido a MRSA são, atualmente, antibióticos glicopeptídicos, como a vancomicina e teicoplanina.


O aparecimento de estirpes resistentes a antibióticos, como é o caso das MRSA e VRSA, e a habilidade de ultrapassar a eficácia dos medicamentos disponíveis no mercado e desenvolver resistências, enfatiza a necessidade urgente de se desenvolverem novas drogas para prevenir e controlar infeções relacionadas com *S. aureus*.


## Procedimento

Recorrendo ao comando Entrez. efetch(), faz-se pesquisas de entradas no Pubmed relacionadas com o organismo de interesse e, utilizando o id do genoma do *S. aureus* no NCBI (BA000018.3), retira-se (sob o formato GenBank) a sequência completa do genoma, bem como anotações contendo informações tais como o organismo e taxonomia, as fontes e referências, contendo também informações acerca dos genes e CDS, tais como os seus nomes, localização na sequência, tradução e chaves a usar nas bases de dados pertinentes ao trabalho. Guarda-se estes dados automaticamente num ficheiro GenBank. As anotações deste ficheiro encontram-se em anexo.


Com o propósito de limitar o número de genes a utilizar no resto do procedimento, de modo a ser proveitoso em termos de custo computacional e de tempo, para além de obtermos já os genes essenciais que será útil para o passo seguinte, recorremos ao modelo do repositório do *S. aureus* N315 e fizemos uma simulação de genes essenciais, com a função objetivo R_biomass_SA_8a, obtendo 168 genes essenciais. Colocam-se os nomes desses genes numa lista no Python, e através de um ciclo, comparando os nomes desses genes à anotação de cada CDS, restringimos a informação a ser utilizada a partir daqui.


Para cada um dos genes selecionados procedeu-se à procura de genes homólogos com o genoma humano, através de BLAST (‘blastn’ por serem sequências de nucleótidos). No BioPython, isto é possível fazer remotamente, através do comando NCBIWWW.qblast, sendo a base de dados usada ‘nr’, e a limitação ao genoma humano é feita pela definição do parâmetro entrez_query='txid9606 [ORGN]' (sendo este código a identificação do genoma humano). Para propósito de diminuir o custo computacional, definiu-se o número máximo de alinhamentos, hitlist_size, como 25 (metade do default).


Definiu-se o threshold de 0.05 para o e-value, em que alinhamentos cujos e-values sejam menores que esse limite são considerados como alinhamentos significativos (portanto, genes com homologia no genoma humano), e valores acima do limite assumem-se como alinhamentos sem homologia. Encontraram-se 101 genes essenciais sem homologia no genoma humano (em anexo).


Para verificar a essencialidade dos genes, verificou-se se estes estão presentes no DEG, que é uma base de dados de genes essenciais. Verificou-se uma grande discrepância entre os genes essenciais obtidos pelo Optflux e os indicados nesta base de dados, mas de modo a garantir análise a genes dos quais se tinha a certeza da sua essencialidade, selecionou-se os genes essenciais presentes tanto no Optflux e DEG, sem homologia nos humanos, sendo o número destes genes apenas 16 (em anexo).


Selecionaram-se 3 genes da lista de genes essenciais descrita no parágrafo anterior: SA1259, SA0997 e SA0457 e as suas respetivas traduções (Dihydrofolate reductase, Glutamate racemase e UDP-N-acetylglucosamine), e para cada destes, retiraram-se as suas informações do Uniprot (tais como a sua sequência de aminoácidos, referências e outras informações pertinentes). Para além disso, também foram lidos os dados do NCBI CDD, em que se obtém os domínios conservados de cada proteína; foram usados o Phobius e Boctopus para encontrar domínios transmembranares alfa e beta, respetivamente; PDB para obter informações estruturais de cada proteína; LocTree3 para a localização sub-celular e CBS para deteção de locais de ligação (por fosforilação); tudo isto para aumentar a compreensão acerca de cada proteína, o que será útil para definição de drogas terapêuticas.


<br />

```python
#bibliotecas utilizadas
from Bio import Entrez #Para aceder ao NCBI
from Bio.Blast import NCBIXML
from Bio.Blast import NCBIWWW
from Bio import ExPASy
from Bio import SwissProt
#%%
#pesquisa automatizada de artigos no Pubmed de interesse
Entrez.email = "...@aaa.org"
handle = Entrez.efetch(db="pubmed", id="28588310,29111265,23519164", retmode="xml")
search_n315=Entrez.read(handle)
 #%%
#genes essenciais presentes na simulação de genes essenciais do modelo do repositório do Optflux
# com função objetivo R_biomass_SA_8a 
genes_essenciais=['SA0923','SA0924','SA0925','SA0926','SA0920','SA0921','SA0922','SA0938','SA0937','SA0916',
                  'SA0917','SA0918','SA0919','SA0912','SA0913','SA0915','SA0910','SA0911','SA1938','SA0842',
                  'SA0843','SA0965','SA1297','SA1299','SA1298','SA1177','SA1052','SA2027','SA1065','SA1397',
                  'SA1150','SA2127','SA1165','SA1164','SA1166','SA2136','SA2186','SA2287','SA1197','SA1199',
                  'SA2288','SA1088','SA1089','SA1571','SA0486','SA0244','SA0487','SA1461','SA1585','SA0375',
                  'SA1228','SA1348','SA1227','SA1229','SA1587','SA0376','SA0134','SA1586','SA1347','SA1226',
                  'SA1589','SA0016','SA1588','SA1104','SA1346','SA1558','SA0347','SA1439','SA2406','SA0344',
                  'SA0345','SA0346','SA0592','SA0472','SA0593','SA0473','SA0594','SA0474','SA1205','SA1202',
                  'SA2412','SA0596','SA1201','SA0597','SA1204','SA1203','SA2413','SA1494','SA2341','SA2465',
                  'SA1496','SA2464','SA1493','SA1250','SA1492','SA2467','SA2466','SA1259','SA2468','SA2347',
                  'SA0176','SA0177','SA2470','SA2471','SA0178','SA0179','SA1352','SA1115','SA1244','SA2333',
                  'SA2456','SA1245','SA1487','SA2334','SA1126','SA1368','SA1858','SA1735','SA1731','SA1860',
                  'SA1749','SA0419','SA1865','SA1861','SA1864','SA1863','SA0506','SA1959','SA0865','SA1608',
                  'SA1729','SA1728','SA0512','SA1965','SA0996','SA0997','SA0756','SA1724','SA0994','SA0995',
                  'SA1651','SA0683','SA1650','SA1652','SA1412','SA0693','SA1309','SA0457','SA0458','SA1427',
                  'SA1669','SA1306','SA1301','SA1422','SA1545','SA1424','SA0549','SA0547','SA0669','SA0548',
                  'SA0670','SA0793','SA0439','SA0794','SA0795','SA1523','SA0796','SA1522']#optflux
genes_essenciais=sorted(genes_essenciais) #ordena alfabeticamente os genes
Entrez.email='lol@lul.com'
handle=Entrez.efetch(db='nucleotide', rettype='gb', retmode='text', id='BA000018.3') #retira o genoma do S.aureus N315
aureus=SeqIO.read(handle, 'genbank') #lê o genoma
SeqIO.write(aureus,"saureusn315.gbk","genbank") #escreve num ficheiro o genoma
handle.close()
#%%
#imprime dados de anotação do genoma
print(aureus.id)
print(aureus.name)
print()
print(aureus.description)
print()
print(aureus.annotations)
print()
print(aureus.dbxrefs)
 #%%
#guarda as regiões codificantes dos genes essenciais do optflux, através da informações dos qualifiers dos features de cada CDS 
pos_genes_essenciais=[]
genesaureus=[]
for i in range(len(aureus.features)):
    if aureus.features[i].type=='CDS':
        essencial=False
        pos=0
        for g_essencial in genes_essenciais:
            if g_essencial in str(aureus.features[i].qualifiers["note"]):
                essencial=True
                pos=genes_essenciais.index(g_essencial)
                break
        if essencial:
            genesaureus.append(i)
            pos_genes_essenciais.append(pos)
genes=[]
for i in genesaureus:
    genes.append(aureus.features[i].extract(aureus.seq))
 #%%
#para cada CDS imprime os qualifires
for i in pos_genes_essenciais:
    print(aureus.features[i])
    print(aureus.features[i].qualifiers)
 #%%
E_VALUE_THRESH = 0.05 # assume-se que o threshold sobre o qual a homologia é insignificativa é 0.05
blast_records_BA=[]
#%%
#para cada gene essencial do Optflux, faz-se o Blast contra o genoma humano (txid9606 [ORGN]) e guarda o alinhamento
for gene in genes:
    result_handle= NCBIWWW.qblast("blastn","nr",gene._data,entrez_query='txid9606 [ORGN]', hitlist_size=25)
    blast_records=NCBIXML.read(result_handle)
    try:
        print ('Homologia.')
        blast_records_BA.append((blast_records.alignments))
    except:
        print('Blast não possível.')
        blast_records_BA.append('Blast não possível.')
#%%      
# verifica-se se os alinhamentos são significativos (com valor de evalue/expect menor que 0.05)
genes_sem_homologia=[]
for i in range(len(blast_records_BA)):
    print(genes_essenciais[i], end= " ")
    try:
       print(blast_records_BA[i][0].hsps[0].expect)
       if blast_records_BA[i][0].hsps[0].expect>E_VALUE_THRESH:
           genes_sem_homologia.append([genes_essenciais[i],blast_records_BA[i][0].hsps[0].expect]) #i+84
    except:
        print("SEM RESULTADOS.")
 #%%
#verifica-se se os genes essenciais obtidos com o Optflux batem certo com os presentes na BD DEG 
genes_DEG=open("genes_essenciais_DEG.txt","r")
genes_essenciais_deg=genes_DEG.readlines()
genes_DEG.close()
g_essenciais_deg=[]
for i in range(len(genes_essenciais_deg)):
    if genes_essenciais_deg[i].find("/")==-1:
        genes_essenciais_deg[i]=genes_essenciais_deg[i].strip('\n')
        g_essenciais_deg.append(genes_essenciais_deg[i])
    else:
        genes_essenciais_deg[i]=genes_essenciais_deg[i].strip('\n')
        temp=genes_essenciais_deg[i].split("/")
        for j in temp:
            g_essenciais_deg.append(j)
genes_essenciais_deg=g_essenciais_deg
 #%%
genes_optflux_gb=[]
for i in range(len(genes_essenciais_deg)):
    for j in range(len(aureus.features)):
        if aureus.features[j].type=="CDS":
            if (str(genes_essenciais_deg[i]) in str(aureus.features[j].qualifiers["note"]) or
                str(genes_essenciais_deg[i]) in str(aureus.features[j].qualifiers["gene"])):
                    genes_optflux_gb.append([genes_essenciais_deg[i],j])
                    break
#%%
genes_deg=[]
for i in genes_optflux_gb:
    ind = i[1]
    string=(str(aureus.features[ind].qualifiers["note"]))
    ind_2=string.find("ORFID:")
    genes_deg.append(string[ind_2+6:ind_2+12])
 #%%
deg_optflux=[]
genes_s_homologia_optflux=[]
#for i in result_s_homologia: genes_s_homologia_optflux.append(i[0])
for i in genes_sem_homologia: genes_s_homologia_optflux.append(i[0])
genes_s_homologia_optflux=sorted(genes_s_homologia_optflux)
for i in genes_deg:
    if i in genes_s_homologia_optflux:deg_optflux.append(i)
    
#%%
#   Genes sem homologia com humanos, presentes no Optflux e DGE
#   Gene    e-value     hiperligação no UniProt
# 'SA0179', 0.309651  https://www.uniprot.org/uniprot/P60296
# 'SA0457', 0.267056 https://www.uniprot.org/uniprot/Q7A7B4
# 'SA0924', 0.14368 https://www.uniprot.org/uniprot/P99162
# 'SA0997', 0.638148 https://www.uniprot.org/uniprot/P63638
# 'SA1104', 0.205002 https://www.uniprot.org/uniprot/Q7A5Y4
# 'SA1177', 1.87107 https://www.uniprot.org/uniprot/P99161
    
#['SA1204', 0.281624 https://www.uniprot.org/uniprot/P66987
# 'SA1259', 0.372101 https://www.uniprot.org/uniprot/P99079
# 'SA1492', 0.889915 https://www.uniprot.org/uniprot/P64334
# 'SA1494', 0.244495 https://www.uniprot.org/uniprot/P64341
# 'SA1728', 0.429948 https://www.uniprot.org/uniprot/P99150
# 'SA2027', 0.588693 https://www.uniprot.org/uniprot/P99062
# 'SA2406'] 0.391652 https://www.uniprot.org/uniprot/A0A0H3JNT1
# 'SA1522', 0.0679558    
# 'SA1346', 0.0967526
 #%%
# lê o ficheiro do Proteome do S. aureus N315
aureus_prot=SeqIO.parse("uniprot-proteome_UP000000751.fasta", "fasta")
#%%
aa_fasta=[]
for i in aureus_prot:
    aa_fasta.append(i)
 #%%
#retira os dados das proteínas do Uniprot
#https://www.uniprot.org/uniprot/P99079
#SA1259  Dihydrofolate reductase  
# annotation score: 3/5
handle = ExPASy.get_sprot_raw('P99079')
record = SwissProt.read(handle)
 #%%
info=[record.accessions, record.annotation_update, record.comments, record.created, record.cross_references, record.data_class, record.description, record.entry_name, record.features, record.gene_name, record.host_organism, record.host_taxonomy_id, record.keywords, record.molecule_type, record.organelle, record.organism, record.organism_classification, record.protein_existence, record.references, record.seqinfo, record.sequence, record.sequence_length, record.sequence_update, record.taxonomy_id]
for i in range(len(info)):
    print(dir(record)[i+26],":" ,end=" ")
    if info[i] != record.references: 
        print(info[i])
        print()
    else:
        for k in range(len(info[i])):
            print()
            print("Authors:" ,end=" ")
            print(info[i][k].authors)
            print("comments:" ,end=" ")
            print(info[i][k].comments)
            print("location:" ,end=" ")
            print(info[i][k].location)
            print("number: ", end=" ")
            print(info[i][k].number)
            print("positions:" ,end=" ")
            print(info[i][k].positions)
            print("references:" ,end=" ")
            print(info[i][k].references)
            print("title:" ,end=" ")
            print(info[i][k].title)
            print()
 #%%
#retira os dados do NCBI CDD para essa proteína
handle=Entrez.efetch(db='protein', rettype='gb', retmode='text', id='P99079')
Dhf_reductase=SeqIO.read(handle, 'genbank')
SeqIO.write(Dhf_reductase,"Dihydrofolate_reductase.gbk","genbank")
handle.close()
for i in Dhf_reductase.features:
    print(i)
    print(i.qualifiers)
     
#%%
#retira os dados das proteínas do Uniprot
# SA0997 Glutamate racemase
# https://www.uniprot.org/uniprot/P63638
# annotation score: 2/5
handle = ExPASy.get_sprot_raw('P63638')
record = SwissProt.read(handle)
 #%%
info=[record.accessions, record.annotation_update, record.comments, record.created, record.cross_references, record.data_class, record.description, record.entry_name, record.features, record.gene_name, record.host_organism, record.host_taxonomy_id, record.keywords, record.molecule_type, record.organelle, record.organism, record.organism_classification, record.protein_existence, record.references, record.seqinfo, record.sequence, record.sequence_length, record.sequence_update, record.taxonomy_id]
for i in range(len(info)):
    print(dir(record)[i+26],":" ,end=" ")
    if info[i] != record.references: 
        print(info[i])
        print()
    else:
        for k in range(len(info[i])):
            print()
            print("Authors:" ,end=" ")
            print(info[i][k].authors)
            print("comments:" ,end=" ")
            print(info[i][k].comments)
            print("location:" ,end=" ")
            print(info[i][k].location)
            print("number: ", end=" ")
            print(info[i][k].number)
            print("positions:" ,end=" ")
            print(info[i][k].positions)
            print("references:" ,end=" ")
            print(info[i][k].references)
            print("title:" ,end=" ")
            print(info[i][k].title)
            print()
            
#%%
#retira os dados do NCBI CDD para essa proteína 
handle=Entrez.efetch(db='protein', rettype='gb', retmode='text', id='P63638')
Glu_racemase=SeqIO.read(handle, 'genbank')
SeqIO.write(Glu_racemase,"Glutamate_racemase.gbk","genbank")
handle.close()
for i in Glu_racemase.features:
    print(i)
    print(i.qualifiers)
    
 #%%
#retira os dados das proteínas do Uniprot
#'SA0457', UDP-N-acetylglucosamine 
#https://www.uniprot.org/uniprot/Q7A7B4
# annotation score: 5/5
 
handle = ExPASy.get_sprot_raw('Q7A7B4')
record = SwissProt.read(handle)
 #%%
info=[record.accessions, record.annotation_update, record.comments, record.created, record.cross_references, record.data_class, record.description, record.entry_name, record.features, record.gene_name, record.host_organism, record.host_taxonomy_id, record.keywords, record.molecule_type, record.organelle, record.organism, record.organism_classification, record.protein_existence, record.references, record.seqinfo, record.sequence, record.sequence_length, record.sequence_update, record.taxonomy_id]
for i in range(len(info)):
    print(dir(record)[i+26],":" ,end=" ")
    if info[i] != record.references: 
        print(info[i])
        print()
    else:
        for k in range(len(info[i])):
            print()
            print("Authors:" ,end=" ")
            print(info[i][k].authors)
            print("comments:" ,end=" ")
            print(info[i][k].comments)
            print("location:" ,end=" ")
            print(info[i][k].location)
            print("number: ", end=" ")
            print(info[i][k].number)
            print("positions:" ,end=" ")
            print(info[i][k].positions)
            print("references:" ,end=" ")
            print(info[i][k].references)
            print("title:" ,end=" ")
            print(info[i][k].title)
            print()
            
 #%%
#retira os dados do NCBI CDD para essa proteína 
handle=Entrez.efetch(db='protein', rettype='gb', retmode='text', id='Q7A7B4')
glmU=SeqIO.read(handle, 'genbank')
SeqIO.write(glmU,"BifunctionalProteinGlmU.gbk","genbank")
handle.close()
for i in glmU.features:
    print(i)
    print(i.qualifiers)
```
<br />

# Proteínas Escolhidas
<br />

## Dihidrofolato redutase


Dihidrofolato redutase (DHFR), com id de acessão P99079 na base de dados Protein do NCBI e no UniProt (com score de anotação 3/5 no UniProt), é proveniente do gene folA, com id de acessão SA1259. Esta enzima é importante no metabolismo do folato, que cataliza a reação representada na figura 1. Esta é uma reação essencial para a síntese de novo de glicina e purina, e para a síntese de precursores de DNA, como a timina. A estrutura tridimensional da enzima encontra-se representada na figura 2, em que se encontra ligado o NADP+.
<br />
<br />
<img src="https://upload.wikimedia.org/wikipedia/commons/thumb/c/ce/DHFR_Reaction_Scheme.png/308px-DHFR_Reaction_Scheme.png">

##### Figura 1 – Reação catalisada pela enzima dihidrofolato redutase, EC 1.5.1.3

<br />

<img src="%23SA1259  Dihydrofolate reductase/DHFR 3D.PNG" width="400">

##### Figura 2 – Imagem tridimensional da dihidrofolato redutase, obtida do PDB através do código de acessão 6E4E. Apesar desta enzima ser proveniente de uma estirpe de *S. aureus* diferente da N315, a sua sequência aminoacídica é a mesma

<br />

Com recurso a biopython, importamos o ficheiro genbank da enzima dihidrofolato redutase. Ainda com o biopython, é possível extrair informações deste ficheiro, nomeadamente das features, que contém informações relevantes sobre regiões da proteína. Neste caso, verificamos que a enzima tem um comprimento de 159 aminoácidos e nas features é possível identificar os locais de ligação dos substratos (NADPH e dihidrofolato) na enzima. Os valores apresentados correspondem aos aminoácidos que participam na ligação dos substratos à enzima. O conhecimento dos locais de ligação e do modo como se ligam os substratos é sempre útil no design de drogas, por permitir encontrar substâncias que apresentem padrões de ligação semelhantes, sendo ainda possível os locais da enzima que servirão como potencias alvos, por serem importantes na ligação ao substrato. Todas as anotações do ficheiro do UniProt estão em anexo.
<br />
De seguida procedemos ao estudo da localização, organização estrutural e modificações pós-tradução desta enzima. Com recurso ao LocTree3, previmos a localização sub-celular da DHFR. O resultado encontra-se na figura 3. A enzima será citoplasmática e tendo em conta a função que desempenha, a localização faz algum sentido. Usámos o Phobius e Boctupus para encontrar regiões α-hélice e β-barril transmembranares, respetivamente. Os resultados encontram-se na figura 4 e verificam a inexistência de domínios transmembranares em ambos os casos, que corrobora o resultado obtido pelo LocTree3, de que a enzima será citoplasmática. Recorrendo ao PDB, usado atrás para a estrutura 3D da enzima, podemos avaliar a existência de α-hélices e folhas β. Na figura 5 encontra-se uma lista de ‘features’ da proteína. Na linha ‘Secstruc’ verificamos a existência de 4 hélices e 10 folhas β, corroborado, também, pela observação da estrutura tridimensional.
<br />

<img src="%23SA1259  Dihydrofolate reductase/loc_sub_celular.PNG" width="700">

##### Figura 3 – Localização subcelular da DHFR prevista pelo LocTree3
<br />

<img src="%23SA1259  Dihydrofolate reductase/Dom_TransMemb_alpha.PNG" width="650"> <img src="%23SA1259  Dihydrofolate reductase/Dom_TransMemb_bet.PNG" width="600">

##### Figura 4 – Previsão de domínios α-hélice transmembranares da DHFR prevista pelo Phobius (A) e de domínios β-barril pelo Boctupus (B)
<br />

<img src="%23SA1259  Dihydrofolate reductase/DHFR SecStruct.PNG" width="800">

##### Figura 5 – Lista de features da DHFR, obtida do PDB, com código de acesso 6E4E. As estruturas α-hélices e folhas β encontram-se representadas a vermelho e bege, respetivamente, na linha Secsstruc
<br />

Procedemos à identificação de possíveis modificações pós-tradução, que alteram as características da proteína, alteram polaridade e tamanho dos aminoácidos, afetando a ligação da proteína ao substrato, a localização subcelular pode ser determinada por estas modificações. Estudamos a fosforilação, nos resíduos de serina, treonina e tyrosina, na DHFR através do NetPhosBac. Os resultados obtidos encontram-se na figura 6, em que T é a treonina, S a serina e Y a tirosina. Tendo representadas as possíveis posições de fosforilação, podemos comparar com as posições de ligação aos substratos, obtidos em biopython, pela lista de features da enzima. Verificamos que existe fosforilação nas posições 36, 40, 79, 136 e 137. Por comparação com as features, verificamos que as fosforilações não afetarão os locais de ligação de substrato, pelo que estas modificações não terão impacto no modo de ligação dos substratos. 
<br />

<img src="%23SA1259  Dihydrofolate reductase/DHFR_locais_fosf.png">

##### Figura 6 – Posições de fosforilação na DHFR previstas pelo NetPhosBac
<br />

Os domínios conservados na sequência também foram determinados, com recurso ao ScanProsite e CDD do NCBI. Os motivos geralmente estão relacionados com determinadas funções biológicas. Assim, a determinação de motivos conservados será útil no design de fármacos, por se ter o conhecimento, por exemplo de locais de ligação de metais, substratos ou outras funções associadas aos motivos. Em ambas as ferramentas é detetado o motivo pertencente à superfamília DHFR (dihidrofolato redutase), característico deste tipo de enzimas. Verificamos que há um “match”, significativo, com quase toda a totalidade da enzima. No caso do CDD ainda é possível verificar as posições de interação dos dois substratos com aminoácidos individuais, que correspondem às mesmas posições presentes nas features obtidas com biopython. Estes locais são potenciais alvos de drogas contra a DHFR.


	Foi realizado um BLAST tendo como query a sequência proteica da DHFR e foram escolhidos 10 organismos, que apresentavam homologia com esta enzima. Com as sequências da DHFR dos 10 organismos e de *S. aureus* realizámos um alinhamento múltiplo, do qual foi depois contruída uma árvore filogenética. A árvore encontra-se na figura 7 e verifica-se um elevado grau de relacionamento entre a DHFR de *S. aureus* com as de Glycine max e Arabidopsis taliana, duas espécies vegetais. Apesar destas duas proteínas terem tamanhos muito superiores à de *S. aureus*, 530 e 408 aminoácidos, respetivamente, demonstram uma zona homologa à DHFR, com uma função associada igual. Relativamente às restantes enzimas bacterianas, todas elas têm um tamanho semelhante à DHFR de *S. aureus*, tendo a mesma função. Estes resultados corroboram os obtidos pelo CDD. A construção de árvores filogenéticas também será útil no estudo de potenciais agentes terapêuticos. Será provável que uma droga que atue numa dada enzima, também seja ativa contra enzimas homólogas de outro organismo. Neste caso, não faz sentido o estudo de drogas para as duas espécies de plantas, devido à maior complexidade do organismo e diferenças significativas relativamente à bactéria. Assim sendo, o descobrimento de drogas como sendo ativas em Escherichia coli, Bacillus subtilis, entre outros será útil, pois estes agentes terapêuticos poderão ser ativos contra a DHFR de *S. aureus*.
	
<br />

<img src="%23SA1259  Dihydrofolate reductase/DHFR_tree.PNG">

##### Figura 7 – Árvore filogenética para enzimas homólogas com a DHFR, obtida do alinhamento múltiplo de organismos selecionados com BLAST
<br />
<br />
Trimethoprim é uma droga que actua como um análogo de pirimidina e que perturba a síntese de folato, essencial para a via de síntese de timidina. A inibição da DHPR faz com que o organismo não sintetise os nucleótidos necessários para a replicação de DNA, actuando como uma bactericida. Trimethoprim liga-se à proteína e inibe a redução de ácido dihidrofólico (DHF) para ácido tetrahidrofólico (THF). Sendo que THF é um precursor essencial na via de síntese de timidina, há deste modo uma inibição da síntese de DNA. Adicionalmente, Trimethoprim apresenta uma muito maior afinidade para a DHPR da bactéria do que para a humana, o que faz desta proteína um bom alvo terapêutico. Trimethoprim actua melhor em conjunto com Sulfamethoxazole, uma droga que inibe outra enzima envolvida na mesma via, a dihidropteroato sintetase. Esta combinação das duas drogas funciona melhor do que o Trimethoprim por si só, pois reduz o desenvolvimento de resistência por parte de *S. aureus* a estas drogas[9].
<br />

<img src="%23SA1259  Dihydrofolate reductase/Trimethoprim.jpg">

##### Figura 8 – Trimethoprim
<br />

<br />

## Glutamato racemase
<br />
A glutamato racemase (GLUR), com id de acessão P63638, na base de dados Protein do NCBI e no UniProt (com score de anotação 2/5 no UniProt), é resultado da transcrição e tradução do gene murI, com id de acessão SA0997. A enzima está envolvida no metabolismo do glutamato, essencial para a biossíntese da parede celular em bactérias, formando D-glutamato a partir de L-glutamato, reação representada na figura 9. D-glutamato é um monómero da camada de peptidoglicano, um componente essencial na estrutura da parede celular em bactérias. A conservação do glutamato racemase e a sua essencialidade no crescimento em procariotas faz desta enzima um bom alvo para a descoberta de potenciais drogas. A estrutura tridimensional da GLUR encontra-se representada na figura 10, em que se encontram ligadas duas moléculas de D-glutamato. Por esta estrutura, é possível perceber que a enzima é um homodímero.
<br />

<img src="%23 SA0997 Glutamate racemase/GLUR reaction.PNG">

##### Figura 9 – Reação catalisada pela enzima glutamato racemase, EC 5.1.1.3
<br />

<img src="%23 SA0997 Glutamate racemase/GLUR 3D.PNG" width="400">

##### Figura 10 – Imagem tridimensional do glutamato racemase, obtida do Swiss Model através do código de acessão P63638.
<br />

De novo, recorrendo ao biopython, importamos o ficheiro genbank da enzima glutamato racemase. Tal como atrás, é possível extrair as features e outras informações relevantes. Neste caso, verificamos que a enzima tem um comprimento de 266 aminoácidos e nas features é possível identificar os locais de ligação do substrato (glutamato) na enzima. Os valores aqui encontrados correspondem aos aminoácidos que participam na ligação do substrato à enzima, conhecimento útil no design de drogas. Todas as anotações do ficheiro do UniProt estão em anexo.
<br />
Procedemos ao estudo da localização, organização estrutural e modificações pós-tradução desta enzima. Com recurso ao LocTree3, previmos a localização sub-celular da GLUR, com o resultado representado na figura 11, sendo a enzima citoplasmática. Usámos o Phobius e Boctupus novamente para encontrar regiões α-hélice e β-barril transmembranares, respetivamente. Os resultados encontram-se na figura 12 e verificam a inexistência de domínios transmembranares em ambos os casos, que corrobora o resultado obtido pelo LocTree3, de que a enzima será citoplasmática. Recorrendo ao PDB, usado atrás para a estrutura 3D da enzima, podemos avaliar a existência de α-hélices e folhas β. Na figura 13 encontra-se uma lista de ‘features’ da proteína. Na linha ‘Secstruc’ verificamos a existência de 11 hélices e 10 folhas β, corroborado, apesar de ser difícil de observar, pela estrutura tridimensional.
<br />

<img src="%23 SA0997 Glutamate racemase/loc_sub_celular.PNG" width="700">

##### Figura 11 – Localização subcelular da GLUR prevista pelo LocTree3
<br />

<img src="%23 SA0997 Glutamate racemase/Dom_TransMemb_alpha.PNG" width="650"> <img src="%23 SA0997 Glutamate racemase/Dom_TransMemb_bet.PNG" width="400">

##### Figura 12 – Previsão de domínios α-hélice transmembranares da GLUR prevista pelo Phobius (A) e de domínios β-barril pelo Boctupus (B)
<br />

<img src="%23 SA0997 Glutamate racemase/GLUR SecStruct.PNG" width="800">

##### Figura 13 – Lista de features da GLUR, obtida do PDB, com código de acesso 2JFQ, de uma strain de S. auresus diferente, mas a sequência é igual. As estruturas α-hélices e folhas β encontram-se representadas a vermelho e bege, respetivamente, na linha Secsstruc
<br />

De seguida identificamos locais de possíveis modificações pós-tradução. Estudamos a fosforilação, nos resíduos de serina, treonina e tyrosina, na GLUR através do NetPhosBac, tal como anteriormente. Os resultados obtidos encontram-se na figura 14, em que T é a treonina, S a serina e Y a tirosina. Tendo representadas as possíveis posições de fosforilação, podemos comparar com as posições de ligação aos substratos, obtidos em biopython, pela lista de features da enzima. Verificamos que existe fosforilação nas posições 88, 122, 168, 175, 185, 206, 207, 243 e 256. Por comparação com as features da GLUR, verificamos que a posição 185 está envolvida na ligação ao glutamato, pelo que esta modificação terá impacto no modo de ligação do substrato.
<br />

<img src="%23 SA0997 Glutamate racemase/GluRace_locais_fosf.PNG">

##### Figura 14 – Posições de fosforilação na GLUR previstas pelo NetPhosBac
<br />

Mais uma vez foram determinados os domínios conservados, com recurso ao ScanProsite e CDD do NCBI. Em ambas as ferramentas é detetado o motivo pertencente à superfamília Asp_Glu_race (Asp/Glu/Hydantoin racemase). Pertencem a esta família racemases relacionadas evolucionariamente, que não necessitam de cofatores para a sua atividade enzimática. Verificamos que há um “match” de duas regiões pequenas, típicas nesta superfamília.


Para a filogenia, foi realizado mais uma vez um BLAST para identificar organismos com homologia na sequência da enzima. Foram escolhidos 9 organismos para fazer alinhamentos múltiplos e construir a árvore filogenética, todos bactérias. As respetivas enzimas têm todas a mesma função (glutamato racemase) e tamanhos muito semelhantes. Estes resultados demonstrativos de homologia entre diferentes espécies ajudam a comprovar os resultados de motivos conservados obtidos com o CDD. A árvore filogenética encontra-se na figura 15. Verifica-se uma maior semelhança entre a GLUR de *S. aureus* com as de *Bacillus subtilis* e *Streptococcus pneumoniae*. De facto, a glutamato racemase é uma enzima com a estrutura e modo de funcionamento muito conservados em diversas bactérias[10]. Assim, o descobrimento de fármacos para a GLUR de espécies bacterianas mais bem estudadas possivelmente levará a estes fármacos também serem ativos contra *S. aureus*.

<br />
<img src="%23 SA0997 Glutamate racemase/GLUR_tree.PNG">

##### Figura 15 - Árvore filogenética para enzimas homólogas com a GLUR, obtida do alinhamento múltiplo de organismos selecionados com BLAST
<br />

Embora a glutamato racemase em *S. aureus* não esteja ainda descrita com alvo de drogas terapêuticas conhecidas, já foram identificados em diferentes espécies de bactérias vários tipos de inibidores desta enzima, incluindo “mechanism and substrate-based inhibitors”, inibidores alostéricos, e “large molecule inhibitors”[10]. Um exemplo de inibidores alostéricos que foram identificados são as pirazolopirimidinedionas, que mostraram inibir a glutamato racemase em Helicobacter pylori[11]. Estes inibidores exibem uma inibição reversível, e a ligação da enzima ao substrato é necessária para que haja inibição, pois o local de ligação do inibidor só fica acessível após a ligação do substrato ao centro ativo devido à deslocação da hélice C-terminal. Embora a actividade destes inibidores tenha sido testada para várias espécies, incluindo *S. aureus*, estes mostraram possuir elevada afinidade apenas para a glutamato racemase de H. pylori[12].

<br />

<img src="%23 SA0997 Glutamate racemase/Pirazolopirimidinediona.jpg">

##### Figura 16 - Pirazolopirimidinediona.

<br />

Outro grupo de inibidores de glutamato racemase são as 8-benzil pteridine-6,7-dionas, com atividade inibitória em bactérias Gram-positivas[13]. Estes inibidores têm origem noutro grupo de inibidores, as 9-benzil purinas, que não demonstraram ter atividade contra a glutamato racemase de *S. aureus*, emboram demonstrassem inibir a de espécies próximas como Enterococcus faecalis e Enterococcus faecium[14]. No entanto, as 8-benzil pteridine-6,7-dionas já demonstraram ter um espectro de acção mais alargado, que inclui *S. aureus*.

<br />

<img src="%23 SA0997 Glutamate racemase/Benzil purina e 8-benzil pteridinediona.jpg">

##### Figura 17 -Benzil purina e 8-benzil pteridinediona.
<br />


A elevada conservação e essencialidade da glutamato racemase, em conjunto com os vários inibidores já identificados, prova que esta enzima é um bom potencial alvo terapêutico, embora se tenha de ter em consideração as pequenas diferenças estruturais e bioquímicas da enzima entre diferentes espécies para o desenvolvimento de novas drogas em *S. aureus*.


<br />

## N-acetilglucosamina-1-fosfato uridiltransferase
<br />
N-acetilglucosamina-1-fosfato uridiltransferase (GLMU), com id de acessão Q7A7B4, na base de dados Protein do NCBI e no UniProt (com score de anotação 5/5 no UniProt), é resultado da transcrição e tradução do gene glmU, com id de acessão SA0457, em Nucleotide. Esta enzima é essencial no metabolismo de aminoaçúcares e pode ser também um alvo terapêutico atrativo. GlmU catalisa a formação de uridine-diphospho-N-acetylglucosamine (UDP-GlcNAc), um precursor importante na biossíntese de peptidoglicano e lipopolissacarídeos tanto em bactérias gram-negativas, como gram-positivas. GLMU tem um papel bifuncional a possuir dois centros ativos funcionalmente autónomos: o centro acetiltransferase e o centro uridiltransferase que residem em dois domínios proteicos distintos. A reação de acetiltransferase ocorre no domínio C-terminal (acetiltransferase) e a reação uridiltransferase ocorre no domínio N-terminal (uridiltransferase), figura 18. Para a estrutura tridimensional da GLMU não foi encontrado nenhum modelo. Por isso, realizamos ‘sequence search’ no PDB e escolhemos a sequência mais similar. Optámos pela enzima com id 4AAW (E-value=2.052E-121; Identities=49%; Positives=67%), do organismo S. pneumoniae, com exatamente a mesma função que a GLMU. Esta enzima está representada na figura 19, em que é evidente a presença de um homotrímero. Na parte superior da figura, N-terminal, encontram-se os locais ativos da reação de uridiltransferase; na parte inferior, C-terminal, encontram-se os locais ativos da atividade acetiltransferase.
<br />
<br />
<img src="%23'SA0457'%2C UDP-N-acetylglucosamine/GLMU reaction.PNG">

##### Figura 18 – Reação catalisada pela enzima GLMU, EC 5.1.1.3
<br />

<img src="%23'SA0457'%2C UDP-N-acetylglucosamine/GLMU 3D.PNG" width="400">

##### Figura 19 – Imagem tridimensional da GLMU, obtida do PDB através do código de acessão 4AAW.
<br />

Novamente com o biopython, importamos o ficheiro genbank da enzima glutamato racemase. Extraímos as features e outras informações relevantes. Todas as anotações do ficheiro do UniProt estão em anexo. Verificámos que a enzima tem um comprimento de 450 aminoácidos e nas features é possível identificar os locais de ligação dos substratos (UDP-GlcNAc, acetil-coA, Mg2+) aos locais da enzima. Estes valores correspondem aos aminoácidos que participam na ligação dos substratos à enzima, útil para o design de drogas.


Novamente procedemos ao estudo da localização, organização estrutural e modificações pós-tradução desta enzima. Com o LocTree3, previmos a localização sub-celular da GLMU. O resultado obtido encontra-se na figura 20, sendo a enzima citoplasmática, resultado também confirmado pelas anotações da página UniProt desta enzima. Tal como atrás, usámos o Phobius e Boctupus para encontrar regiões α-hélice e β-barril transmembranares, respetivamente. Os resultados encontram-se na figura 21 e verificam que não existem domínios transmembranares em ambos os casos, que confirma o resultado obtido no LocTree3. Recorrendo ao PDB, podemos avaliar a existência de α-hélices e folhas β, recorrendo à enzima usada atrás para avaliar a estrutura 3D. Na figura 22 encontra-se uma lista de ‘features’ desta enzima. Na linha ‘Secstruc’ verificamos a existência de 15 hélices α e cerca de 41 folhas β. 
<br />

<img src="%23'SA0457'%2C UDP-N-acetylglucosamine/loc_celular.PNG" width="700">

##### Figura 20 – Localização subcelular da GLMU prevista pelo LocTree3
<br />

<img src="%23'SA0457'%2C UDP-N-acetylglucosamine/Dom_TransMemb_alpha.PNG" width="650"> <img src="%23'SA0457'%2C UDP-N-acetylglucosamine/Dom_TransMemb_bet.PNG" width="400">

##### Figura 21 – Previsão de domínios α-hélice transmembranares da GLMU prevista pelo Phobius (A) e de domínios β-barril pelo Boctupus (B)
<br />

<img src="%23'SA0457'%2C UDP-N-acetylglucosamine/GLMU SecStruct.PNG" width="800">

##### Figura 22 – Lista de features da GLMU, obtida do PDB, com código de acesso 4AAW, de S. pneumonieae. As estruturas α-hélices e folhas β encontram-se representadas a vermelho e bege, respetivamente, na linha Secsstruc
<br />

Com recurso ao NetPhosBac estudámos os locais de fosforilação nos resíduos de serina, treonina e tirosina, na GLMU. Os resultados obtidos encontram-se na figura 23, em que T é a treonina, S a serina e Y a tirosina. Comparámos as possíveis posições de fosforilação, com as posições de ligação aos substratos, obtidos em biopython, pela lista de features da GLMU. Verificámos que poderá existir fosforilação nas posições 18, 96, 146, 169, 284, 312, 318, 364, 400 e 434. Por comparação com as features da GLMU, verificamos nenhuma destas posições corresponde a locais de ligação a substrato ou cofatores, pelo que as fosforilações não terão impacto na interação entre enzima e substrato.
<br />

<img src="%23'SA0457'%2C UDP-N-acetylglucosamine/GLMU_locais_fosf.PNG">

##### Figura 23 – Posições de fosforilação na GLMU previstas pelo NetPhosBac
<br />

Utilizámos o ScanProsite e CDD do NCBI, para determinar os motivos conservados na GLMU. No caso do ScanProsite, foram detetadas duas regiões conservadas, ‘LuxR-type’ HTH e ‘Hexapeptide-repeat containing-transferases’. No primeiro caso foi identificado um domínio de ligação ao DNA. No entanto, tendo em conta as funções já descritas desta enzima e o nível de confiança obtido, este motivo não estará associado à enzima. A segunda região conservada está associada à atividade de transferase. Tendo em conta a descrição atrás da atividade de acetiltransferase da GLMU, este motivo estará associado a essa atividade, visto a sua posição ser também no C-terminal. No CDD foi encontrada uma superfamília GlmU. Verificámos que há uma homologia forte com as enzimas pertencentes a esta família, todas elas com as mesmas atividades.


A GLMU foi ainda analizada por BLAST, de modo a encontrar organismos que partilham homologia. Foram escolhidos 10 organismos, todos eles bactérias. Com a sequência das enzimas dos organismos foram realizados alinhamento múltiplo e foi contruída a árvore filogenética, representada na figura 24. As enzimas dos 10 organismos apresentaram todos sequências e funções semelhantes à GLMU de *S. aureus*, resultados que reforçam as funções associadas à enzima. Da árvore verificamos ainda que há um relacionamento forte entre *S. aureus* e *Bacillus subtilis*, também uma bactéria gram-positiva, e ainda com as enzimas de *Clostridioides difficile* e *Streptococcus pneumoniae*, também gram-positivas. Com este conhecimento, será provável que uma droga desenvolvida para um destes organismos seja ativa contra a enzima de *S. aureus*, devido à elevada homologia, pelo que estudos nestes organismos poderão ser interessantes.



<img src="%23'SA0457'%2C UDP-N-acetylglucosamine/GLMU_tree.PNG">

##### Figura 24 – Árvore filogenética para enzimas homólogas com a GLMU, obtida do alinhamento múltiplo de organismos selecionados com BLAST



Foi identificada uma molécula sintética de pequeno tamanho que inibe a GlmU em Haemophilus influenzae, que ocupa um centro alostérico adjacente ao local de ligação do substrato GlcNAc-1-P, previnindo, deste modo, a ocorrência de rearranjos estruturais necessários para que haja a reação catalisada pela enzima[15]. Não se verificou, no entanto, atividade contra a GlmU de bactérias Gram-positivas, incluindo a de *S. aureus*. A GlmU destes organismos contém uma substituição num resíduo no centro alostérico relativamente à GlmU de Haemophilus influenzae, o que impede a interação entre o inibidor e a enzima. Contudo, estes resultados sugerem que o centro alostérico pode ser utilizado para criar compostos com uma melhor afinidade para ortólogos de GlmU, abrindo assim as portas para o desenvolvimento de uma nova classe de drogas em *S. aureus* tendo esta enzima com alvo terapêutico.

<br />

<img src="%23'SA0457'%2C UDP-N-acetylglucosamine/Molécula inibidora da GlmU.jpg">

##### Figura 25 - Molécula inibidora da GlmU.

<br />
#################################(REGULAÇÃO)#################################

# Referências

1.	 Masalha M, Borovok I, Schreiber R, Aharonowitz Y, Cohen G (2001) Analysis of transcription of the *Staphylococcus aureus* aerobic class Ib and anaerobic class III ribonucleotide reductase genes in response to oxygen. Journal of Bacteriology. 183 (24): 7260–72. doi:10.1128/JB.183.24.7260-7272.2001. PMC 95576. PMID 11717286.


2.	 Chambers HF (2001) The changing epidemiology of *Staphylococcus aureus*? Emerging Infectious Diseases. 7 (2): 178–82. doi:10.3201/eid0702.010204. PMC 2631711. PMID 11294701.


3.	 Jevons MP (1961) Celbenin-resistant staphylococci. BMJ. 1 (5219): 124–5. doi:10.1136/bmj.1.5219.124-a.


4.	 Deurenberg, R. H., & Stobberingh, E. E. (2008) The evolution of *Staphylococcus aureus*. Infection, genetics and evolution, 8(6), 747-763.


5.	 Khan, M. F. (2017) Brief History of *Staphylococcus aureus*: A Focus to Antibiotic Resistance. EC Microbiology, 5, 36-39.


6.	 Blot SI, Vandewoude KH, Hoste EA, Colardyn FA (2002) Outcome and attributable mortality in critically III patients with bacteremia involving methicillin-susceptible and methicillin-resistant *Staphylococcus aureus*. Archives of Internal Medicine. 162 (19): 2229–35. doi:10.1001/archinte.162.19.2229. PMID 12390067.


7.	 Hiramatsu K, Hanaki H, Ino T, Yabuta K, Oguri T, Tenover FC (1997) Methicillin-resistant *Staphylococcus aureus* clinical strain with reduced vancomycin susceptibility. The Journal of Antimicrobial Chemotherapy. 40 (1): 135–6. doi:10.1093/jac/40.1.135. PMID 9249217.


8.	 Carter AP, Clemons WM, Brodersen DE, Morgan-Warren RJ, Wimberly BT, Ramakrishnan V (2000) Functional insights from the structure of the 30S ribosomal subunit and its interactions with antibiotics. Nature. 407 (6802): 340–8. doi:10.1038/35030019. PMID 11014183.


9.	 Drugbank: Trimethoprim. https://www.drugbank.ca/drugs/DB00440


10.	Fisher SL (2008) Glutamate racemase as a target for drug discovery. Microbial Biotechnology. 1 (5): 345–60. doi:10.1111/j.1751-7915.2008.00031.x


11.	Lundqvist T, Fisher SL, Kern G, Folmer RHA, Xue Y, Newton DT, et al. (2007) Exploitation of structural and regulatory diversity in glutamate racemases. Nature. 447: 817–22


12.	de Jonge BLM, Kutschke A, Uria-Nickelsen M, Kamp HD, Mills SD (2009) Pyrazolopyrimidinediones are selective agents for Helicobacter pylori that suppress growth through inhibition of glutamate racemase (MurI). American Society for Microbiology. 53 (8): 3331–6. doi:10.1128/AAC.00226-09.


13.	Breault G, Eyermann CJ, Comita-Prevoir J, Geng B, Petrichko R. (2007) Hit to lead studies: exploring 8-benzyl pteridine 6,7-diones as inhibitors of glutamate racemase (MurI) in Gram positive bacteria. 47th ICAAC Meeting, Chicago, IL, USA, Abs: F1-336.


14.	Geng B, Breault G, Comita-Prevoir J, Petrichko R, Eyermann JC, Doig P, Gorseth E, Noonan B (2008) Bioorg. Med. Chem. Lett., 18, 4368.


15.	Mochalkin I, Lightle S, Narasimhan L, Bornemeier D, Melnick M, Vanderroest S, McDowell L (2008) Structure of a small-molecule inhibitor complexed with GlmU from Haemophilus influenzae reveals an allosteric binding site. Protein Science. 17 (3): 577–82. doi: 10.1110/ps.073271408.


