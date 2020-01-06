# Reproduce table 1
set.seed(47)
for(gamma in seq(0,1,length.out = 5)){
  for(n in c(200,400)){
    reptable(gamma = gamma, n = n,rep = 1000,clusterN = 100 )
  }
}

# Reproduce table 2
for(gamma in seq(0,1,length.out = 3)){
  for(n in c(400)){
    reptable(gamma = gamma, n = n,rep = 1000,clusterN = 100 )
  }
}
