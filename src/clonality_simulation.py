#! /usr/bin/python

############################################################################
## A agent-based model to simulate 3D tumor growth and metastasis         ##
## Copyright: Zheng Hu, Stanford University                               ##
############################################################################

import sys,os,math,random
import numpy as np
from collections import Counter
import sets
from virtualTumor3D import *


############main script#################
repl = int(sys.argv[1])       #replicate of simulations
seeding_size = int(sys.argv[2])  #monoclonal seeding: seeding_size=1; polyclonal seeding: seeding_size=10; 

rd = 60 #radius of the 3D space
mut_rate = 0.6 #per cell-division mutation rate
p_adv_rate = pow(10,-5) #rate of advantageous mutations in primary tumor
m_adv_rate = pow(10,-5) #rate of advantageous mutations in metastasis
p_scoef = pow(10,random.uniform(-3,-1)) #selective coefficient of advantageous mutations in primary tumor
m_scoef = pow(10,random.uniform(-3,-1)) #selective coefficient of advantageous mutations in metastasis

p_birth_rate = random.uniform(0.55,0.65) #birth probability for each cell at each generation in primary tumor
m_birth_rate = random.uniform(0.55,0.65) #birth probability for each cell at each generation in metastasis

td = int(2*pow(10,random.uniform(0,4)))  #timing of dissemination
npub = 50
tissue_purity = 0.6 # tumor purity
seq_depth = 100 # mean sequencing depth
deme_size = 5000 # deme size
final_tumor_size = pow(10,9)
final_deme_number = final_tumor_size/deme_size


###################################
flist = ["0"]
mut_lineage = ['0']          # the lineage tracer
mut_index = 0

###Neutral growth of primary tumor
#neu_cells,mut_index,mut_lineage = initiateFirstDeme_neutral(deme_size,mut_lineage,flist,mut_index,p_birth_rate)  #the growth of the fisrt deme from single transformed cell

###Selective growth of primary tumor
neu_cells,adv_cells,mut_index,mut_lineage = initiateFirstDemeSel_v1(deme_size,mut_lineage,flist,mut_index,p_scoef,p_adv_rate,p_birth_rate)  #the growth of the fisrt deme from single transformed cell

################
space = createLattice(rd)
space[(rd,rd,rd)].present = 1
space[(rd,rd,rd)].neutral = list(neu_cells)
space[(rd,rd,rd)].advant = list(adv_cells)
current_keys = [(rd,rd,rd)]
surface_keys = [(rd,rd,rd)]
current_deme_number =1 #current tumor size as # of demes
surface_deme_number =1 #current tumor size as # of demes
generation = 0

while current_deme_number < final_deme_number:
    new_keys = []
    for w in range(0,surface_deme_number):
        skey = random.choice(surface_keys)
        if space[skey].present == 1:
            nei_sites = neighbor26(skey)
            empty_neis = [key for key in nei_sites if space[key].present == 0]                    # empty neighbor sites
            
            if len(empty_neis) > 0:
                rand_prob = random.random()
                if rand_prob < 1-math.exp(-len(empty_neis)*0.25):   #when the deme will be chosen for growth and fission
                    pre_neu = list(space[skey].neutral)
                    pre_adv = list(space[skey].advant)
                    post_neu_1,post_neu_2,post_adv_1,post_adv_2,mut_index,mut_lineage = demeGrowthFissionSel_v1(pre_neu,pre_adv,mut_lineage,mut_index,current_deme_number,p_scoef,p_adv_rate,p_birth_rate)
                    newkey = random.choice(empty_neis)
                    new_keys += [newkey]
                    space[skey].neutral = list(post_neu_1)
                    space[newkey].neutral = list(post_neu_2)
                    space[skey].advant = list(post_adv_1)
                    space[newkey].advant = list(post_adv_2)
                    space[newkey].present = 1
                    current_keys += [newkey]
                    current_deme_number += 1

                    if current_deme_number == td:
                        if seeding_size == 1:
                            anc_list = list(post_neu_2+post_adv_2)
                            met_progenitor = random.sample(anc_list,seeding_size)
                        else:
                            current_surface_keys = list(surface_keys+new_keys)
                            all_cells = all_mutation_ids(space,current_surface_keys)
                            met_progenitor = random.sample(all_cells,seeding_size)
                        
    
    ###update surface
    surface_update = list(surface_keys+new_keys)
    surface_keys = []
    for fkey in surface_update:
        neisites = neighbor26(fkey)
        random.shuffle(neisites)
        for key in neisites:
            if space[key].present == 0:
                surface_keys += [fkey]
                break
    surface_deme_number = len(surface_keys)
    generation += 1
    

#################
pquad1,pquad2,pquad3,pquad4 = [],[],[],[]
for pky in surface_keys:
    if pky[0] > rd and pky[1] > rd and pky[2] > rd:
        pquad1 += [pky]
    if pky[0] < rd and pky[1] < rd and pky[2] < rd:
        pquad2 += [pky]
    if pky[0] < rd and pky[1] > rd and pky[2] > rd:
        pquad3 += [pky]
    if pky[0] > rd and pky[1] < rd and pky[2] < rd:
        pquad4 += [pky]


sample_p1 = random.choice(pquad1)
sample_p2 = random.choice(pquad2)
sample_p3 = random.choice(pquad3)
sample_p4 = random.choice(pquad4)
p4samples = [sample_p1,sample_p2,sample_p3,sample_p4]

ptissue1 = bulkTissueSampling(space,p4samples[0],5)
ptissue2 = bulkTissueSampling(space,p4samples[1],5)
ptissue3 = bulkTissueSampling(space,p4samples[2],5)
ptissue4 = bulkTissueSampling(space,p4samples[3],5)

ppool4 = ptissue1+ptissue2+ptissue3+ptissue4

pvaf1 = seqProcessing(space,ptissue1,mut_lineage,2,seq_depth,tissue_purity)
pvaf2 = seqProcessing(space,ptissue2,mut_lineage,2,seq_depth,tissue_purity)
pvaf3 = seqProcessing(space,ptissue3,mut_lineage,2,seq_depth,tissue_purity)
pvaf4 = seqProcessing(space,ptissue4,mut_lineage,2,seq_depth,tissue_purity)

pvaf_pool4 = seqProcessing(space,ppool4,mut_lineage,2,seq_depth*4,tissue_purity)

#####################################
####growth of met###
####################################

#met_mut_lineage = list(mut_lineage)
#met_mut_index = int(mut_index)

#m_neu_cells,mut_index,mut_lineage = initiateFirstDeme_neutral(deme_size,mut_lineage,met_progenitor,mut_index,m_birth_rate)  #the growth of the fisrt deme from single transformed cell
m_neu_cells,m_adv_cells,mut_index,mut_lineage = initiateFirstDemeSel_v1(deme_size,mut_lineage,met_progenitor,mut_index,m_scoef,m_adv_rate,m_birth_rate)  #the growth of the fisrt deme from single transformed cell

#print m_adv_cells
#print met_progenitor

mspace = createLattice(rd)
mspace[(rd,rd,rd)].present = 1
mspace[(rd,rd,rd)].neutral = list(m_neu_cells)
mspace[(rd,rd,rd)].advant = list(m_adv_cells)

met_current_keys = [(rd,rd,rd)]
met_surface_keys = [(rd,rd,rd)]
met_current_deme_number =1                                                 #current number of demes
met_surface_deme_number =1                                                 #current number of demes

        
while met_current_deme_number < final_deme_number:
    met_new_keys = []
    for w in range(0,met_surface_deme_number):
        skey = random.choice(met_surface_keys)
        if mspace[skey].present == 1:
            nei_sites = neighbor26(skey)
            empty_neis = [key for key in nei_sites if mspace[key].present == 0]                    # empty neighbor sites
            
            if len(empty_neis) > 0:
                rand_prob = random.random()
                if rand_prob < 1-math.exp(-len(empty_neis)*0.25):   #when the deme will be chosen for growth and fission
                    pre_neu = list(mspace[skey].neutral)
                    pre_adv = list(mspace[skey].advant)
                    post_neu_1,post_neu_2,post_adv_1,post_adv_2,mut_index,mut_lineage = demeGrowthFissionSel_v1(pre_neu,pre_adv,mut_lineage,mut_index,met_current_deme_number,m_scoef,m_adv_rate,m_birth_rate)
                    newkey = random.choice(empty_neis)
                    met_new_keys += [newkey]
                    mspace[skey].neutral = list(post_neu_1)
                    mspace[newkey].neutral = list(post_neu_2)
                    mspace[skey].advant = list(post_adv_1)
                    mspace[newkey].advant = list(post_adv_2)
                    mspace[newkey].present = 1
                    met_current_keys += [newkey]
                    met_current_deme_number += 1
    
    ###update surface
    met_surface_update = list(met_surface_keys + met_new_keys)
    met_surface_keys = []
    for fkey in met_surface_update:
        neisites = neighbor26(fkey)
        random.shuffle(neisites)
        for key in neisites:
            if mspace[key].present == 0:
                met_surface_keys += [fkey]
                break
    met_surface_deme_number = len(met_surface_keys)
        

#################
mquad1,mquad2,mquad3,mquad4 = [],[],[],[]
for pky in met_surface_keys:
    if pky[0] > rd and pky[1] > rd and pky[2] > rd:
        mquad1 += [pky]
    if pky[0] < rd and pky[1] < rd and pky[2] < rd:
        mquad2 += [pky]
    if pky[0] < rd and pky[1] > rd and pky[2] > rd:
        mquad3 += [pky]
    if pky[0] > rd and pky[1] < rd and pky[2] < rd:
        mquad4 += [pky]

sample_m1 = random.choice(mquad1)
sample_m2 = random.choice(mquad2)
sample_m3 = random.choice(mquad3)
sample_m4 = random.choice(mquad4)
m4samples = [sample_m1,sample_m2,sample_m3,sample_m4]

mtissue1 = bulkTissueSampling(mspace,m4samples[0],5)
mtissue2 = bulkTissueSampling(mspace,m4samples[1],5)
mtissue3 = bulkTissueSampling(mspace,m4samples[2],5)
mtissue4 = bulkTissueSampling(mspace,m4samples[3],5)

mpool4 = mtissue1+mtissue2+mtissue3+mtissue4

mvaf1 = seqProcessing(mspace,mtissue1,mut_lineage,2,seq_depth,tissue_purity)
mvaf2 = seqProcessing(mspace,mtissue2,mut_lineage,2,seq_depth,tissue_purity)
mvaf3 = seqProcessing(mspace,mtissue3,mut_lineage,2,seq_depth,tissue_purity)
mvaf4 = seqProcessing(mspace,mtissue4,mut_lineage,2,seq_depth,tissue_purity)

mvaf_pool4 = seqProcessing(mspace,mpool4,mut_lineage,2,seq_depth*4,tissue_purity)

MAF_file = open("Tumor3D_MRS_8samples_c"+str(seeding_size)+"_"+str(repl)+".txt","w")
MAF_file.write("mut_id"+" "+"public"+" "+"pdepth1"+" "+"pvaf1"+" "+"pdepth2"+" "+"pvaf2"+" "+"pdepth3"+" "+"pvaf3"+" "+"pdepth4"+" "+"pvaf4"+" "+"mdepth1"+" "+"mvaf1"+" "+"mdepth2"+" "+"mvaf2"+" "+"mdepth3"+" "+"mvaf3"+" "+"mdepth4"+" "+"mvaf4")
MAF_file.write("\n")

MAF_pool4 = open("Tumor3d_MRS_pooled_c"+str(seeding_size)+"_"+str(repl)+".txt","w")
MAF_pool4.write("mut_id"+" "+"public"+" "+"pdepth"+" "+"pvaf"+" "+"mdepth"+" "+"mvaf")
MAF_pool4.write("\n")


for k in range(0,npub):
    pdepth,pvaf = pubMutGenerator(8,2,seq_depth,tissue_purity)
    MAF_file.write("0"+" "+"1"+" "+str(pdepth[0])+" "+str(pvaf[0])+" "+str(pdepth[1])+" "+str(pvaf[1])+" "+str(pdepth[2])+" "+str(pvaf[2])+" "+str(pdepth[3])+" "+str(pvaf[3])+" "+str(pdepth[4])+" "+str(pvaf[4])+" "+str(pdepth[5])+" "+str(pvaf[5])+" "+str(pdepth[6])+" "+str(pvaf[6])+" "+str(pdepth[7])+" "+str(pvaf[7]))
    MAF_file.write("\n")

    pdepth,pvaf = pubMutGenerator(2,2,seq_depth*4,tissue_purity)
    MAF_pool4.write("0"+" "+"1"+" "+str(pdepth[0])+" "+str(pvaf[0])+" "+str(pdepth[1])+" "+str(pvaf[1]))
    MAF_pool4.write("\n")


muts_all = sets.Set(pvaf1.keys()) | sets.Set(pvaf2.keys()) | sets.Set(pvaf3.keys()) | sets.Set(pvaf4.keys()) |sets.Set(mvaf1.keys()) | sets.Set(mvaf2.keys()) |sets.Set(mvaf3.keys()) | sets.Set(mvaf4.keys())
muts_all_pool4 = sets.Set(pvaf_pool4.keys()) | sets.Set(mvaf_pool4.keys())

pabsent1 = muts_all-sets.Set(pvaf1.keys())
pabsent2 = muts_all-sets.Set(pvaf2.keys())
pabsent3 = muts_all-sets.Set(pvaf3.keys())
pabsent4 = muts_all-sets.Set(pvaf4.keys())
mabsent1 = muts_all-sets.Set(mvaf1.keys())
mabsent2 = muts_all-sets.Set(mvaf2.keys())
mabsent3 = muts_all-sets.Set(mvaf3.keys())
mabsent4 = muts_all-sets.Set(mvaf4.keys())

pabsent_pool4 = muts_all_pool4 - sets.Set(pvaf_pool4.keys())
mabsent_pool4 = muts_all_pool4 - sets.Set(mvaf_pool4.keys())

pvaf1 = missingDepth(pvaf1,pabsent1,seq_depth)
pvaf2 = missingDepth(pvaf2,pabsent2,seq_depth)
pvaf3 = missingDepth(pvaf3,pabsent3,seq_depth)
pvaf4 = missingDepth(pvaf4,pabsent4,seq_depth)
mvaf1 = missingDepth(mvaf1,mabsent1,seq_depth)
mvaf2 = missingDepth(mvaf2,mabsent2,seq_depth)
mvaf3 = missingDepth(mvaf3,mabsent3,seq_depth)
mvaf4 = missingDepth(mvaf4,mabsent4,seq_depth)

pvaf_pool4 = missingDepth(pvaf_pool4,pabsent_pool4,seq_depth*4)
mvaf_pool4 = missingDepth(mvaf_pool4,mabsent_pool4,seq_depth*4)

for mt in list(muts_all):
    MAF_file.write(str(mt)+" "+"0"+" "+str(pvaf1[mt][0])+" "+str(pvaf1[mt][1])+" "+str(pvaf2[mt][0])+" "+str(pvaf2[mt][1])+" "+str(pvaf3[mt][0])+" "+str(pvaf3[mt][1])+" "+str(pvaf4[mt][0])+" "+str(pvaf4[mt][1])+" "+str(mvaf1[mt][0])+" "+str(mvaf1[mt][1])+" "+str(mvaf2[mt][0])+" "+str(mvaf2[mt][1])+" "+str(mvaf3[mt][0])+" "+str(mvaf3[mt][1])+" "+str(mvaf4[mt][0])+" "+str(mvaf4[mt][1]))
    MAF_file.write("\n")

for mt in list(muts_all_pool4):
    MAF_pool4.write(str(mt)+" "+"0"+" "+str(pvaf_pool4[mt][0])+" "+str(pvaf_pool4[mt][1])+" "+str(mvaf_pool4[mt][0])+" "+str(mvaf_pool4[mt][1]))
    MAF_pool4.write("\n")

