load("data.codon_nbs")
load("iterate_vals_base")
source('/Users/keren/Dropbox/mtDNA_tree_Build_17/debug_functions2_base_thesis.R')

get_nb_model = function(codon_position,codon){
  tmp = get_data(iterate_vals,j = codon_position, q = codon)
  curr_data = tmp[[1]]; 
  if (is.na(curr_data)==FALSE){
    tmp = prepare_data_rm_cols(curr_data)
    curr_data = tmp[[1]]; curr_model_names = tmp[[2]]
    model.nb = find_model.nb.data(curr_data,curr_model_names)
    # model.p = find_model.P.data(curr_data, curr_model_names)
  }else{
    model.nb = NA
  }
  return(model.nb)
}

# The winning model partitions according to codon position and codon, takes the Right and Left 
# neighbours and the direction of replication as explaining variables and ignores all other 
# variables. Hence it is composed of 3*64=192 sub-models. We provide below default values for the
# variables i,k,p,w,t,u,m,n and a, all determeined by the winning model. Values for j and q vary
# (j = 1:3, q = 1:64) to describe all the sub-models that compose the winning model.

i = 5; # Output options:
# 1:syn = synonymous substitutions, 2:non_syn = non synonymous substitutions,
# 3:transitions, 4:transversions, 5:y = all substitutions
k = 3 # Directionality: 1:all genes except for ND6, 2: ND6 (ND6 has opposite directionality)
# 3: directionality is an explaining variable, 4: directionality is not included in the regression.
p = 5 # CG position (1-3, not CG, first of CG pair, second of CG pair), 
# 4:CG position is an explaining variable, 5:CG position is not included in the regression.
w = 23 # Amino acids (1-21), 22: Amino acids are an explaining variable, 
# 23: Amino acids are not included in the regression.
t = 5; u=5 # Right (t) and Left (u) neighbours (1-4), 5:explaining variable, 6:not included
m = 15 # Genes (m) (1-13 coding genes), 14:explaining variable, 15:not included
n = 5 # Aggregated genes (m) (1-3), 4:explaining variable, 5:not included
a = 6 # Input nucleotide (1-4), 5:explaining variable, 6:not included

#j = codon position (1-3), 4:codon position is an explaining variable, 
# 5:codon position is not included in the regression.
#q = codon (1-64), 65:explaining variable, 66:not included.

model.nb = list()
count = 0
for (j in 1:3){
  for(q in 1:64){
    count = count + 1
    print(count)
    model.nb[[count]] = get_nb_model(codon_position = j,codon = q)
  }
}

