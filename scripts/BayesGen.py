
from scipy import stats
import numpy.random

#let us generate N_REG regulators, from normal distribution mixtures with means 0,1 and std set to param STD

#Noise parameters, reasonable are from (0,1) interval
STD=0.002
PROB_ERR=0.0001

# number of regulators
N_REG=2
#number of target genes
N_TARG=4

# number of observations for each regulator  combination
N_OBS=2



noised=[]
bins=[]
for i in range(2**N_REG):
    bin_v=bin(i)[2:]
    bin_padded="0"*(N_REG-len(bin_v))+bin_v
    for j in range(N_OBS):
        bins.append([int(x) for x in bin_padded])
        noised.append([numpy.random.randn()*STD+int(x) for x in bin_padded])

#noisy and
def noisy_and(*args,**kwds):
    #print args[0]
    if sum(args[0])==len(args[0]): #AND==1
        if numpy.random.rand()<PROB_ERR: #error is made
            return 0
        else: # no error
            return 1
    else: # AND==0
        if numpy.random.rand()<PROB_ERR: #error is made
            return 1
        else: # no error
            return 0

#generating noisy and targets
#print len(bins)

targ_res=[]
for i in range(2**N_REG*N_OBS): #observation number
    res_i=[]
    for j in range(N_TARG): #target gene number
        res_i.append(numpy.random.randn()*STD+int(noisy_and(bins[i][:j])))#target j takes first j regulators
        #print bins[i][:j], noisy_and(bins[i][:j]), res_i[-1]
    targ_res.append(res_i)

#print the results
print "#regulators", " ".join(["gene_%s"%i for i in range(N_REG)])
hdr="gene\t"+"\t".join(["exp_%d"%i for i in range(2**N_REG*N_OBS)])
print hdr

for i in range(N_REG):
    print "gene_%d\t"%i + "\t".join(["%g"%noised[j][i] for j in range(2**N_REG*N_OBS)])
for i in range(N_TARG):
    print "target_gene_%d\t"%i + "\t".join(["%g"%targ_res[j][i] for j in range(2**N_REG*N_OBS)])
        
