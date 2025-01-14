library(edgeR)
library(pryr) # per utilizzare unenclose() -> funzioni non utilizzato l'environment generale


#correlation with method Pearson
corr <- unenclose(function(dat){
  cor(t(dat),method = "pearson")
})

# Fisher Transformation formula
ft <- unenclose(function(x){
  return (1/2 * log((1+x)/(1-x)))
})

# Z-scores formula
z <- unenclose(function(x,y,n1,n2){
  v1 = 1/(n1-3)
  v2 = 1/(n2-3)
  zf = x - y
  vf = v1 + v2
  
  return (zf/sqrt(vf))
})

# Questa funzione non viene rinchiusa in un enclose() perchè fa utilizzo di 
# funzioni esterne, tuttavia non utilizza alcun parametro esterno.
evaluate_zscores <- function(matrix.condition1,matrix.condition2){
  
  # similarity matrix condition 1
  sim.condition1 <- corr(matrix.condition1)
  
  # similarity matrix condition 2
  sim.condition2 <- corr(matrix.condition2)
  
  # remove diagonal
  diag(sim.condition1) <- diag(sim.condition1) <- 0
  diag(sim.condition2) <- diag(sim.condition2) <- 0
  
  
  # apply function to the respective matrices
  FT.sim.condition1 = apply(sim.condition1, 2, ft)
  FT.sim.condition2 = apply(sim.condition2, 2, ft)
  
  
  n1 <- ncol(matrix.condition1)
  n2 <- ncol(matrix.condition2)
  
  # apply function to the respective tumor types
  zscores <- mapply(z,as.data.frame(FT.sim.condition1),
                    as.data.frame(FT.sim.condition2),
                    n1,
                    n2
  )
  
  rownames(zscores) <- colnames(zscores)
  
  return(zscores)
}


differential_coexpression_network <-  unenclose(function(z,threshold){
  
  library(network)
  # Calcolare la adjacency.matrix basandosi sulla threshold fornita in input
  ## Se z >= t assegnare 1, sennò andare all'altro ifelse
  ## l'altro ifelse fa: se z <= -t assegnare 1 sennò 0
  adjacency.matrix <- ifelse(z > (threshold), 1, ifelse(z < (-threshold), -1, 0))  

  #----------------------------------------------------------------------------#
  
  # Calcolare il degree per ogni gene
  degree <- rowSums(abs(adjacency.matrix))
  
  ## Creare un dataframe che fornisca il nome (enseble) per ogni gene ed il suo degree
  degree <- as.data.frame(cbind(colnames(adjacency.matrix),degree))
  
  ## Rinominare le colonne del dataframe per maggiore chiarezza
  colnames(degree) <- c("gene","degree")
  
  ## Trasformare la colonna degree da stringa in numeri interi
  degree$degree <- as.integer(degree$degree)
  
  #----------------------------------------------------------------------------#
  # Genes with a degree below 1 were excluded 
  degree.zero <- degree[degree$degree > 0,]
  #----------------------------------------------------------------------------#
  
  # Escludere dalla adjacency.matrix i geni con degree pari a 0.
  
  ## Trovare il nome (enseble) dei geni che hanno degree pari a 0
  excluded.genes <- setdiff(degree$gene, degree.zero$gene)
  
  if (length(excluded.genes) != 0){
    ## Trovare l'indice (enseble) dei geni che hanno degree pari a 0
    idx <- which(colnames(adjacency.matrix) %in% excluded.genes)
    
    ## Rimuovere i geni  con degree pari a 0 sfruttando il loro  indice 
    adjacency.matrix<- adjacency.matrix[-idx,-idx]
  }
  
  
  

  #----------------------------------------------------------------------------#
  # Costruire il network utilizzando l' adjacency matrix
  net <- network::network(adjacency.matrix, matrix.type="adjacency",ignore.eval = T, directed = F)
  
  
  
  #----------------------------------------------------------------------------#
  #Definire l'output
  output <- list( adjacency.matrix = adjacency.matrix,
                  network = net,
                  geni.esclusi = excluded.genes
                )
  
  
  return(output)
})