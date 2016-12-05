library('huge')
library('MASS')

# Set directory
setwd('./desktop/stat241_proj/mullets')

# Load data for data w/ LD  
n_sim = 100
num_locs = 15
num_gen = 100 
p = num_gen * num_locs
allele_trajs = matrix(nrow=n_sim, ncol=p)
for (i in 1:n_sim){
  sim_i = i - 1 
  df = read.csv(paste('./data/neutral_withLD/sim.s0.iter',as.character(sim_i), '.csv', sep=''), sep='\t', header=TRUE)
  df = df[,-1] # First column is just ids 
  allele_trajs[i, ] = as.vector(as.matrix(df))
}

# Do graphical model selection 
trans_data = huge.npn(allele_trajs)
mb_select = huge(trans_data, method='glasso', sym='or')
adj_mat = mb_select$path[10][[1]]

# Let's just look at marginals for now 
tm0 = c(1, 201, 401, 601, 801) + 5
tm1 = tm0 + 1
tm2 = tm1 + 1
tm3 = tm2 + 1
tm4 = tm3 + 1
tm5 = tm4 + 1
#all_times = c(tm0, tm1, tm2, tm3, tm4, tm5)
all_times = c(tm0, tm1, tm2)
subset_adj = adj_mat[all_times, all_times]
subset_adj = as.data.frame(as.matrix(subset_adj))

names0 = c('A0', 'B0', 'C0', 'D0', 'E0')
names1 = c('A1', 'B1', 'C1', 'D1', 'E1')
names2 = c('A2', 'B2', 'C2', 'D2', 'E2')
#names3 = c('A3', 'B3', 'C3', 'D3', 'E3')
#names4 = c('A4', 'B4', 'C4', 'D4', 'E4')
#names5 = c('A5', 'B5', 'C5', 'D5', 'E5')
#all_names = c(names0, names1, names2, names3, names4, names5)
all_names = c(names0, names1, names2)
rownames(subset_adj) = all_names
colnames(subset_adj) = all_names

# W/ selection
n_sim = 100
num_locs = 15
num_gen = 100 
p = num_gen * num_locs
allele_trajs = matrix(nrow=n_sim, ncol=p)
for (i in 1:n_sim){
  sim_i = i - 1 
  df = read.csv(paste('./data/selection_withLD/sim.s100.iter',as.character(sim_i), '.csv', sep=''), sep='\t', header=TRUE)
  df = df[,-1] # First column is just ids 
  allele_trajs[i, ] = as.vector(as.matrix(df))
}

trans_data = huge.npn(allele_trajs)
mb_select = huge(trans_data, method='mb', sym='or')
adj_mat = mb_select$path[10][[1]]

# time t = 5
tm0 = c(1, 201, 401, 601, 801) + 5
tm1 = tm0 + 1
tm2 = tm1 + 1
tm3 = tm2 + 1
all_times = c(tm0, tm1, tm2)
subset_adj = adj_mat[all_times, all_times]
subset_adj = as.data.frame(as.matrix(subset_adj))

all_names = c(names0, names1, names2)
rownames(subset_adj) = all_names
colnames(subset_adj) = all_names

# time t = 49
tm0 = c(1, 201, 401, 601, 801) + 49 
tm1 = tm0 + 1
tm2 = tm1 + 1
all_times = c(tm0, tm1, tm2)
subset_adj = adj_mat[all_times, all_times]
subset_adj = as.data.frame(as.matrix(subset_adj))

all_names = c(names0, names1, names2)
rownames(subset_adj) = all_names
colnames(subset_adj) = all_names
