#! /usr/bin/python

############################################################################
## A agent-based model to simulate 3D tumor growth and metastasis         ##
## Copyright: Zheng Hu, Stanford University                               ##
############################################################################

import sys,os,math,random
import numpy as np
from collections import Counter
import sets

class deme():
    def __init__(self):
        self.present= 0     # whether the deme is empty or occupied: 0-empty;1-occupied
        self.neutral = []   # the background cells during tumor growth
        self.advant = []    # the cells with advantageous driver mutations


def createLattice(d):
    """
    A function to create 3D cubic lattice with length of 2d where each site contains a empty deme
    """
    lattice = {}
    for x in range(0,2*d+1):
        for y in range(0,2*d+1):
            for z in range(0,2*d+1):
                lattice[(x,y,z)] = deme()
    return lattice


def neighbor26((a,b,c)):
    """
    Moore neighbourhood: 26 neighbour sites of (a,b,c)
    """
    neighbor = [(a+i,b+j,c+k)
                for i in [-1,0,1]
                for j in [-1,0,1]
                for k in [-1,0,1]
                if not (i==0 and j==0 and k==0)]
    
    return neighbor


def neighbor6((a,b,c)):
    """
    von Neumann neighbourhood: 6 neighbour sites of (a,b,c)
    """
    neighbor = [(a-1, b, c),(a+1, b, c),(a, b-1, c),(a, b+1, c),(a, b, c-1),(a, b, c+1)]
    return neighbor


def localNeighbor((a,b,c),r):
    """
    A function to search the local neighbour sites of (a,b,c) within an area of radius r in the 3D cubic lattice
    """
    neighbor = []
    for x in range(-r,r+1):
        for y in range(-r,r+1):
            for z in range(-r,r+1):
                if pow(x,2)+pow(y,2)+pow(z,2) < pow(r+1,2):
                    neighbor += [(a+x,b+y,c+z)]
    return neighbor


def traceLineage(mlineage,mutid):
    """
    A function to obtain the mutational lineage of a cell from the mutation id of the most recently occurred mutation in the cell. 
    For example, the input id (most recently occurred mutation) of target cell is "100" and the output is "1-12-35-56-100", which is the mutation lineage of the cell
    
    mlineage - the list that could help recover the mutational lineage given the most recent mutation id of a lineage
    mutid - the mutation ID of the most recently occurred mutation in the cell
    """
    recent_muts = mutid.split(',') # it is possible that multiple mutations occur during in a cell division. For instance, the mutation id of most recently occurred mutations is "100,101"
    recent_muts = [int(t) for t in recent_muts]
    first_mut = recent_muts[0] # the first mutation id in a multi-mutation event
    trace = []
    while first_mut > 0:
        trace += recent_muts
        recent_muts = mlineage[first_mut].split(',')
        #print recent_muts
        recent_muts = [int(t) for t in recent_muts]
        first_mut = recent_muts[0]
    return trace

    
def lowerORupper(value):
    lower_int = int(value)
    upper_int = lower_int+1
    if random.random() < value-lower_int:
        return upper_int
    else:
        return lower_int


def initiateFirstDeme_neutral(maxsize,lineage,initial_id,current_id,birth_prob):
    """
    The growth of the initial deme via a neutral process

    maxsize - size limit of a deme
    lineage - a list that stores the lineage information of mutations
    """
    neu_list = [str(x) for x in initial_id]
    n1 = len(neu_list)
    current_deme_size = n1
    while current_deme_size < maxsize:
        neu_divcells =  int(n1*birth_prob+1) #number of dividing cells in this generation
        neu_list = random.sample(neu_list,neu_divcells)*2
        n1 = len(neu_list)
        current_deme_size = n1

        new_mut1 = np.random.poisson(mut_rate*n1)
        mut_assig1 = Counter(np.random.choice(n1,new_mut1))
        for x1 in mut_assig1.keys():
            nmut = mut_assig1[x1]
            new_mut1 = range(current_id+1,current_id+1+nmut)
            mut_str = ",".join(map(str,new_mut1))
            #if nmut > 1:
            #    for t in new_mut1:
            #        multi_events[str(t)] = mut_str
            for xn in range(0,nmut):
                current_id += 1
                lineage += [neu_list[x1]]
            neu_list[x1] = mut_str
    
    return neu_list,current_id,lineage


def initiateFirstDemeSel_v1(maxsize,lineage,initial_id,current_id,sfit,advrate,birth_prob):
    """
    The growth of the initial deme via a selection process.

    maxsize - size limit of a deme
    lineage - a list that stores the lineage information of mutations
    """
    neu_list = [str(x) for x in initial_id]
    adv_list = []
    current_deme_size = len(neu_list)
    while current_deme_size < maxsize:
        n1,n2 = len(neu_list),len(adv_list)
        neu_divcells =  int(n1*birth_prob+1) #number of dividing cells in this generation
        neu_list = random.sample(neu_list,neu_divcells)*2
        if n2 > 0:
            adv_divcells = lowerORupper(n2*birth_prob*(1+sfit)) #number of dividing cells in this generation        
            adv_list = random.sample(adv_list,adv_divcells)*2
        
        n1,n2 = len(neu_list),len(adv_list)
        current_deme_size = n1+n2
        if n1 > 0:
            new_mut1 = np.random.poisson(mut_rate*n1)
            mut_assig1 = Counter(np.random.choice(n1,new_mut1))
            for x1 in mut_assig1.keys():
                nmut = mut_assig1[x1]
                new_mut1 = range(current_id+1,current_id+1+nmut)
                mut_str = ",".join(map(str,new_mut1))
                #if nmut > 1:
                #    for t in new_mut1:
                #        multi_events[str(t)] = mut_str
                for xn in range(0,nmut):
                    current_id += 1
                    lineage += [neu_list[x1]]
                neu_list[x1] = mut_str
        if n2 > 0:
            new_mut2 = np.random.poisson(mut_rate*n2)
            mut_assig2 = Counter(np.random.choice(n2,new_mut2))
            for x2 in mut_assig2.keys():
                nmut = mut_assig2[x2]
                new_mut2 = range(current_id+1,current_id+1+nmut)
                mut_str = ",".join(map(str,new_mut2))
                #if nmut > 1:
                #    for t in new_mut2:
                #        multi_events[str(t)] = mut_str
                for xn in range(0,nmut):
                    current_id += 1
                    lineage += [adv_list[x2]]
                adv_list[x2] = mut_str
        
        if random.random() < advrate*n1:
            current_id += 1
            current_n1 = len(neu_list)
            lineage += [str(neu_list[current_n1-1])]
            adv_list += [str(current_id)]
            neu_list = neu_list[0:current_n1-1]
        
    
    return neu_list,adv_list,current_id,lineage



def demeGrowthFissionSel_v1(neu_list,adv_list,lineage,current_id,current_size,sfit,advrate,birth_prob):
    """
    A function to simulate deme growth and fission and keep track of the mutational lineages
    
    assumptions:
    (1) the growth of neutral and advantageous clones are independent
    """
    current_deme_size = len(neu_list)+len(adv_list)
    while current_deme_size < 2*deme_size:
        n1,n2 = len(neu_list),len(adv_list)
        neu_divcells =  lowerORupper(n1*birth_prob) #number of dividing cells in this generation
        neu_list = random.sample(neu_list,neu_divcells)*2
        if n2 > 0:
            adv_divcells =  lowerORupper(n2*birth_prob*(1+sfit)) #number of dividing cells in this generation
            adv_list = random.sample(adv_list,adv_divcells)*2
        
        n1,n2 = len(neu_list),len(adv_list)
        current_deme_size = n1+n2
        if current_size < 5*pow(10,7)/deme_size: # stop mutation occurence when the tumor size is larger than 10^4*5000 = 5*10^7
            if n1 > 0:
                new_mut1 = np.random.poisson(mut_rate*n1)
                mut_assig1 = Counter(np.random.choice(n1,new_mut1))
                for x1 in mut_assig1.keys():
                    nmut = mut_assig1[x1]
                    new_mut1 = range(current_id+1,current_id+1+nmut)
                    mut_str = ",".join(map(str,new_mut1))
                    #if nmut > 1:
                    #    for t in new_mut1:
                    #        multi_events[str(t)] = mut_str
                    for xn in range(0,nmut):
                        current_id += 1
                        lineage += [neu_list[x1]]
                    neu_list[x1] = mut_str
            if n2 > 0:
                new_mut2 = np.random.poisson(mut_rate*n2)
                mut_assig2 = Counter(np.random.choice(n2,new_mut2))
                for x2 in mut_assig2.keys():
                    nmut = mut_assig2[x2]
                    new_mut2 = range(current_id+1,current_id+1+nmut)
                    mut_str = ",".join(map(str,new_mut2))
                    #if nmut > 1:
                    #    for t in new_mut2:
                    #        multi_events[str(t)] = mut_str
                    for xn in range(0,nmut):
                        current_id += 1
                        lineage += [adv_list[x2]]
                    adv_list[x2] = mut_str
            
            if random.random() < advrate*n1:
                current_id += 1
                current_n1 = len(neu_list)
                lineage += [str(neu_list[current_n1-1])]
                adv_list += [str(current_id)]
                neu_list = neu_list[0:current_n1-1]
            
    
    random.shuffle(neu_list)
    if len(neu_list) > 0:
        offspring_neu = np.random.binomial(len(neu_list),0.5)
    else:
        offspring_neu = 0
    neu_list1=neu_list[0:offspring_neu]
    neu_list2=neu_list[offspring_neu:len(neu_list)]
    
    random.shuffle(adv_list)
    if len(adv_list) > 0:
        offspring_adv = np.random.binomial(len(adv_list),0.5)
    else:
        offspring_adv = 0
    adv_list1=adv_list[0:offspring_adv]
    adv_list2=adv_list[offspring_adv:len(adv_list)]
    
    return neu_list1,neu_list2,adv_list1,adv_list2,current_id,lineage


def seqProcessing(sp,sample_keys,mlineage,size_par,mean_depth,purity):
    """
    Model the random sampling process in next-generation sequencing and report the sequencing allele frequencies in a sample of cells
    
    sp- the lattice space constituting the population of demes
    sample_keys- the sampled demes
    size_par- variance parameter in negative-binomial distribution
    mean_depth- the mean depth of the sequencing
    purity- tumor purity
    """
    all_cur_id = []
    all_mut_id = []
    for key in sample_keys:
        smuts = list(sp[key].neutral + sp[key].advant)
        all_cur_id += smuts
    sample_size = 10000 ## the number of cells for sequencing analysis
    sample_id = random.sample(all_cur_id,sample_size)
    id_count = Counter(sample_id)
    for x in id_count.keys():
        xlineage = traceLineage(mlineage,x)
        all_mut_id += xlineage*id_count[x]
    mut_count = Counter(all_mut_id)
    prob_par=size_par*1.0/(size_par+mean_depth)
    sampleAF = {}
    for x in mut_count.keys():
        true_af = mut_count[x]*0.5*purity/sample_size
        if true_af > 0.005:
            site_depth = np.random.negative_binomial(size_par,prob_par)
            if site_depth >= 10:
                var_reads = np.random.binomial(site_depth,true_af)
                seq_af = var_reads*1.0/site_depth
                if var_reads >= 3:
                    sampleAF[str(x)] = (site_depth,seq_af)
                    #sampleAF[str(x)] = seq_af
    return sampleAF


def highMutCCF(sp,sample_keys,mlineage,cutoff):
    """
    The mutations with cff larger than ccf cutoff
    """
    all_cur_id = []
    all_mut_id = []
    for key in sample_keys:
        smuts = list(sp[key].neutral + sp[key].advant)
        all_cur_id += smuts
    sample_size = 10000 ## the number of cells for sequencing analysis
    sample_id = random.sample(all_cur_id,sample_size)
    id_count = Counter(sample_id)
    for x in id_count.keys():
        xlineage = traceLineage(mlineage,x)
        all_mut_id += xlineage*id_count[x]
    mut_count = Counter(all_mut_id)
    sampleCCF = {}
    for x in mut_count.keys():
        true_ccf = mut_count[x]*1.0/sample_size
        if true_ccf > cutoff:
            sampleCCF[str(x)] = true_ccf
    return sampleCCF


def highMutDeme(sp,position,mlineage,cutoff):
    """
    Obtain the high-frequency mutations (vaf>cutoff) in the deme at site "position"
    """
    all_cur_id = sp[position].neutral + sp[position].advant
    #all_cur_id = sp[position].neutral + sp[position].advant +  sp[position].advant2
    all_mut_id = []
    sample_size = 100
    sample_id = random.sample(all_cur_id,sample_size)
    id_count = Counter(sample_id)
    for y in id_count.keys():
        xlineage = traceLineage(mlineage,y)
        all_mut_id += xlineage*id_count[y]
    mut_count = Counter(all_mut_id)
    highAF_muts = []
    for x in mut_count.keys():
        allele_freq = mut_count[x]*1.0/sample_size
        if allele_freq > cutoff:
            highAF_muts += [x]
    return highAF_muts


def pubMutGenerator(n,size_par,mean_depth,purity):
    """
    A function to generate the public mutations (number is "n") and frequencies
    
    n- number of samples
    size_par- variation parameter in the negative binomial distribution
    mean_death- mean seq depth
    """
    prob_par=size_par*1.0/(size_par+mean_depth)
    mean_af = 0.5*purity
    depth_pub = []
    vaf_pub = []
    for k in range(0,n):
        correct = 0
        while correct == 0:
            site_depth = np.random.negative_binomial(size_par,prob_par)
            if site_depth >= 15:
                correct =1
        var_reads = np.random.binomial(site_depth,mean_af)
        site_vaf = var_reads*1.0/site_depth
        depth_pub += [site_depth]
        vaf_pub += [site_vaf]
    return depth_pub,vaf_pub


def locationSampling(region,sample_number,cutoff):
    """
    A function to sampling the locations where the bulk tissue will be sampled locally.
    """
    success = 0
    while success == 0:
        locations = random.sample(region,sample_number)
        repeat = sample_number*(sample_number-1)
        minall = 999
        for x in range(0,repeat):
            rs = random.sample(locations,2)
            min_distance = min([abs(rs[0][0]-rs[1][0]),abs(rs[0][1]-rs[1][1]),abs(rs[0][2]-rs[1][2])])
            if min_distance < minall:
                minall = min_distance
        if min_distance > 2*cutoff:
            success = 1
    return locations


def bulkTissueSampling(sp,location,radius):
    """
    Check if the sampled sites have non-empty demes
    """
    local_region = localNeighbor(location,radius)
    bulk_tissue = []
    for x in local_region:
        if sp[x].present == 1:
            bulk_tissue += [x]
    return bulk_tissue


def lineageDashLink(mlist):
    """
    Transform the mutation lineage from list (e.g [1,3,10,20]) to dash-linked string (e.g. 1-3-10-20)
    """
    if len(mlist) > 0:
        dstring = str(mlist[0])
        for x in mlist[1:len(mlist)]:
            dstring += "-"
            dstring += str(x)
        return dstring
    else:
        return "0"
        
def missingDepth(vafdata,absent_muts,mean_depth):
    """
    Randomly generate the sequencing depth for the mutation-absent sites the across samples
    """
    for x in absent_muts:
        done = 0
        while done == 0:
            missing_depth = np.random.negative_binomial(2,2.0/(2+mean_depth))
            if missing_depth >= 15:
                done = 1
        vafdata[str(x)] = (missing_depth,0)
    return vafdata


def ccfTransform(vaf_dict):
    """
    Transform the allele frequency to cancer cell fraction (ccf)
    """
    ccf_estimation = {}
    for x in vaf_dict.keys():
        if vaf_dict[x][1] > 0.5:
            ccf_estimation[x] = 1
        else:
            ccf_estimation[x] = vaf_dict[x][1]*2
    return ccf_estimation


def statPrivate(p_mut,m_mut,cutoff):
    p_private = sets.Set(p_mut.keys())-sets.Set(m_mut.keys())
    m_private = sets.Set(m_mut.keys())-sets.Set(p_mut.keys())
    p_stat,m_stat = 0,0
    for x in p_private:
        if p_mut[x] >= cutoff:
            p_stat += 1
    for y in m_private:
        if m_mut[y] >= cutoff:
            m_stat += 1
    return p_stat,m_stat


def all_mutation_ids(sp,sample_keys):
    allkeys = []
    for x in sample_keys:
        allkeys += sp[x].neutral
        allkeys += sp[x].advant
    
    return allkeys

