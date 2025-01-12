# Prende in input la matrice x e le lib.size (che siano sate calcolate o passate in input)
calcFactorQuantile.NEW <- function (data, lib.size, p = 0.75) 
{
  
  # Creare un vettore di tanti 1 quante le colonne di x/data, ovvero per il numero di sample
  # Ovvero inizializza il vettore dei futuri pesi "f"
  f <- rep_len(1, ncol(data))
  
  # Iterare sulle colonne/sample di "data", ovvero le counts già processate dalla funzione
  # seq_len(ncol(data)) restituisce vettore di interi da 1 alla lunghezza passata in input.
  # In questo caso ncol(data). La successione di numeri è per step di 1. equivalente a seq(1,ncol(data))
  for (j in seq_len(ncol(data))) 
    
    # Riassegnare i pesi f con il valore del 75% percentile del sample (j)
    f[j] <- quantile(data[, j], probs = p)
  
  # Se ci sono dei valori di "f" pari a zero sollevare un warning
  if (min(f) == 0) 
    warning("One or more quantiles are zero")
  
  # Restituire, In R, l'uso esplicito di return() per restituire un valore da una funzione è opzionale.
  # Quando una funzione non utilizza return(), l'output della funzione sarà semplicemente 
  # l'ultimo valore calcolato (o l'ultimo valore "valutato") nel corpo della funzione.
  f/lib.size
}

# Funzione per calcolare i fattori di normalizzazione
calcFactorTMM.NEW <- function (obs, ref, libsize.obs = NULL, libsize.ref = NULL, 
                               logratioTrim = 0.3, sumTrim = 0.05, doWeighting = TRUE, 
                               Acutoff = -1e+10) 
{
  
  # Prendere solo i valori numerici del vettore con le  counts per il sample
  obs <- as.numeric(obs)
  
  # Prendere soli i valori numerici del vettore con le  counts per il sample di riferimento
  ref <- as.numeric(ref)
  
  # Se le libsizes del sample non sono specificate ricalcolarle. ed assegnarle ad n0
  if (is.null(libsize.obs)) 
    
    nO <- sum(obs)
  else nO <- libsize.obs # Se invece sono specificate utilizzare quelle specificate
  
  # Se le libsizes del sample di riferimento non sono specificate ricalcolarle ed assegnarle ad nR
  if (is.null(libsize.ref)) 
    nR <- sum(ref)
  else nR <- libsize.ref  #Se invece sono specificate utilizzare quelle specificate
  
  # Continuare con i calcoli 
  
  # Dipede solo da iesimo sample e sample di riferimento
  logR <- log2((obs/nO)/(ref/nR))
  
  # Dipede solo da iesimo sample e sample di riferimento
  absE <- (log2(obs/nO) + log2(ref/nR))/2
  
  # Dipede solo da iesimo sample e sample di riferimento
  v <- (nO - obs)/nO/obs + (nR - ref)/nR/ref
  
  
  fin <- is.finite(logR) & is.finite(absE) & (absE > Acutoff)
  logR <- logR[fin]
  absE <- absE[fin]
  v <- v[fin]
  if (max(abs(logR)) < 1e-06) 
    return(1)
  
  
  n <- length(logR)
  loL <- floor(n * logratioTrim) + 1
  hiL <- n + 1 - loL
  loS <- floor(n * sumTrim) + 1
  hiS <- n + 1 - loS
  keep <- (rank(logR) >= loL & rank(logR) <= hiL) & (rank(absE) >= 
                                                       loS & rank(absE) <= hiS)
  if (doWeighting) 
    f <- sum(logR[keep]/v[keep], na.rm = TRUE)/sum(1/v[keep], 
                                                   na.rm = TRUE)
  else f <- mean(logR[keep], na.rm = TRUE)
  if (is.na(f)) 
    f <- 0
  2^f
}

TMMbyHand <- function(object,lib.size = NULL,
                      method = "TMM",
                      refColumn = NULL,
                      logratioTrim = 0.3, 
                      sumTrim = 0.05,
                      doWeighting = TRUE,
                      Acutoff = -1e+10,
                      p = 0.75){
  
  # Trasformare dataframe delle conte, nome di riferimento OBJECT, in input in una matrice x
  x <- as.matrix(object)
  
  # La funzione "anyNA(x)" restituisce TRUE se ci sono degli NA nella matrice.
  # In caso ci siano NA nella matrice l'esecuzione della funzione viene interrota  
  if (anyNA(x)) 
    stop("NA counts not permitted")
  
  # Calcolare numero di sample della matrice x. 
  # Siccome i samples vengono specificati come colonne ricavare il numero di samples tramite "ncol()"
  nsamples <- ncol(x)
  
  # Se non sono state calcolate ancora le lib.size, calcolarle dato che poi serviranno
  if (is.null(lib.size)) {
    lib.size <- colSums(x) # Le libsize di un sample è la somma dell'espressione di tutti i suoi geni.
    # Perciò calcolare la libsize di ogni sample tramite colSums(x)
    # che restituisce un vettore contenente il totale per ogni colonna (nsamples X 1)  
  }
  
  # Se le lib.size sono già state calcolate allora fare diversi controlli:
  else {
    
    # Per prima cosa controllare che non ci siano NA nelle libsize. 
    # Se ci sono NA l'esecuzione della funzione viene interrota.  
    if (anyNA(lib.size)) 
      stop("NA lib.sizes not permitted")
    
    # Poi controllare che ci siano tante libsize quanti samples
    # Se le libsize sono in numero diverso dal numero di samples:
    if (length(lib.size) != nsamples) {
      
      # Se c'è più di una libsize sollevare un warning
      if (length(lib.size) > 1L) 
        warning("normLibSizes: length(lib.size) doesn't match number of samples", 
                call. = FALSE)
      
      # Risolvere il problema con le libsize, rendendo le libsizes della lunghezza giusta.
      # La funzione rep_len(x=c(1,2),times=3) estende il vettore x per farlo arrivare alla lunghezza times
      # Fa questo ripetendo in ordine gli elementi rep_len(x=c(1,2),times=6) sarà uguale a c(1,2,1,2,1,2)
      lib.size <- rep_len(x=lib.size, times=nsamples)
    }
  }
  

  
  # Booleano, TRUE se una riga (gene) ha tutti valori pari a zero. False altrimenti
  allzero <- .rowSums(x > 0, nrow(x), nsamples) == 0L 
  # il . prima di .rowSums è una convenzione per indicare che .rowSums all'interno di questa funzione
  # è una funzione privata che non deve essere utilizzata dall'utente
  # 0L è una convenzione per indicare che si sta lavorando con gli INTERI
  
  # Se c'è almeno un TRUE in allzero, quindi se c'è almeno un gene con tutti valori pari a zero
  # eliminare quel o quei geni dalla matrice x (escluderli dal processo)
  if (any(allzero)) 
    x <- x[!allzero, , drop = FALSE] # Prende solo quelli con almeno un valore diverso da zero
    
  # Se c'è solo un sample (colonna) o non ci sono geni (righe) nella matrice x method = None
  # con un matrice di queste dimensioni infatti non è possibile applicare le procedure come TMM
  if (nrow(x) == 0 || nsamples == 1) 
    method = "none"
  
  
  # Creare uno switch, ovvero a seconda del metodo specificato intraprendere l'azione
  # Switch eliminato in questa funzione perchè voglio usare solo il TMM
  

  # Se non è stata specificata nessuna refColumn calcolarsela
  # tramite la funzione
  if (is.null(refColumn)) {
    
    # Calcolare il valore di f per ogni sample.
    # Questo valore per ogni sample è pari: 
    # 1. per ogni sample prendere tutte le espressioni del trascrittoma dei geni
    # 2. prendere il valore del 75% percentile
    # 3. Dividere il valore del 75% per la lib.size del gene.
    # 4. Fare questo per ogni sample
    
    f75 <- suppressWarnings(calcFactorQuantile.NEW(data = x, 
                                                   lib.size = lib.size, p = 0.75))
    
    # Se la mediana di tutti i valori f75 (uno per ogni sample) è inferiore a 1e-20.
    # Selezionare il sample di riferimento tramite il criterio: which.max(colSums(sqrt(x)))
    if (median(f75) < 1e-20) {
      # which.max(colSums(sqrt(x))) Si compone di diversi passaggi:
      # 1. "sqrt(X)" # fare la radice quadrata di ogni valore della matrice
      # 2. "colSums(sqrt(X))" ottenere un vettore che contenga le somme per ogni colonna/sample 
      # 3. "which.max(colSums(sqrt(x)))" prendere l'indice del sample 
      # che ha il valore massimo tra le somme di ogni colonna/sample
      
      refColumn <- which.max(colSums(sqrt(x)))
    } else { # Se invece mediana di tutti i valori f75 (uno per ogni sample) è maggiore o uguale di 1e-20.
      # Selezionare la colonna/sample di riferimento tramite il criterio:  which.min(abs(f75 - mean(f75)))
      
      # 1. "f75 - mean(f75)" ad ogni valore f75 sottrarre la media dei valori f75
      # 2. "abs(f75 - mean(f75))" fare il valore assoluto del vettore risultante per ogni elemento
      # 3. "which.min(abs(f75 - mean(f75)))" prendere l'indice del sample che ha il 
      # valore minimo tra i valori del vettore risultante del punto 2.
      refColumn <- which.min(abs(f75 - mean(f75)))
    }
  }
  
  
  # Inizializzare un vettore f, tanti NA quanti il numero di samples
  f <- rep_len(NA_real_, nsamples)
  
  # Calcolare i valore di f per ogni sample.
  # Viene calcolato richiamando la funzione calcFactorTMM
  for (i in 1:nsamples)
    f[i] <- calcFactorTMM.NEW(obs = x[,i], # Selezionare l'iesimo sample/colonna
                              ref = x[, refColumn], # Passare la colonna/sample di riferimento
                              libsize.obs = lib.size[i],  # Passare la lib size dell'iesimo sample 
                              libsize.ref = lib.size[refColumn], # Passare la libsize del sample di riferimento
                              logratioTrim = logratioTrim,  # Parametro
                              sumTrim = sumTrim, # Parametro
                              doWeighting = doWeighting,# Parametro
                              Acutoff = Acutoff # Parametro
                              )
  
  divisore.normalizzazione <- exp(mean(log(f)))
  
  f <- f/exp(mean(log(f)))
  names(f) <- colnames(x)
  

  output <- list(
    counts = object,
    
    fattori.normalizzazione = f,
    divisore = divisore.normalizzazione,
    
    refColumn = refColumn, 
    sample.riferimento = x[, refColumn,drop=FALSE],
    libsize.sample.riferimento = lib.size[refColumn],
    lib.size = lib.size,
    
    dim = dim(x),
    filtered.counts = x,
    filtro = allzero
  )

  return(output)
}

TMMbyHand.forAnotherDataset <- function(object,
                              refColumn,
                              divisore.normalizzazione,
                              lib.size.refColumn=NULL
                              ){
  logratioTrim = 0.3
  sumTrim = 0.05
  doWeighting = TRUE
  Acutoff = -1e+10
  p = 0.75
  

  # geni comuni per confronti tra dataset con geni diversi
  common.genes <- intersect(rownames(object),
                            rownames(refColumn)
                            )
  
  object <- object[common.genes,,drop=FALSE]
  refColumn <- refColumn[common.genes,,drop=FALSE]

  lib.size.refColumn <- sum(refColumn)
  
  # Trasformare dataframe delle conte, nome di riferimento OBJECT, in input in una matrice x
  x <- as.matrix(object)
  
  # La funzione "anyNA(x)" restituisce TRUE se ci sono degli NA nella matrice.
  # In caso ci siano NA nella matrice l'esecuzione della funzione viene interrota  
  if (anyNA(x)) 
    stop("NA counts not permitted")
  
  # Calcolare numero di sample della matrice x. 
  # Siccome i samples vengono specificati come colonne ricavare il numero di samples tramite "ncol()"
  nsamples <- ncol(x)
  
  # Le libsize di un sample è la somma dell'espressione di tutti i suoi geni.
  lib.size <- colSums(x) 


  # Booleano, TRUE se una riga (gene) ha tutti valori pari a zero. False altrimenti
  #
  allzero <- .rowSums(x > 0, nrow(x), nsamples) == 0L 
  # il . prima di .rowSums è una convenzione per indicare che .rowSums all'interno di questa funzione
  # è una funzione privata che non deve essere utilizzata dall'utente
  # 0L è una convenzione per indicare che si sta lavorando con gli INTERI
  
  # Se c'è almeno un TRUE in allzero, quindi se c'è almeno un gene con tutti valori pari a zero
  # eliminare quel o quei geni dalla matrice x (escluderli dal processo)
  if (any(allzero)) 
    x <- x[!allzero, , drop = FALSE] # Prende solo quelli con almeno un valore diverso da zero
  
    # I geni tra x/counts e refColumn devono essere gli stessi, quindi prendere i geni comuni tra entrambe
    # le lib sizes non devono essere aggiornate. (TMM è costruita cosi, guardare il file R originale per conferma)
    refColumn <- refColumn[rownames(x),,drop=FALSE]
  
  # Se c'è solo un sample (colonna) o non ci sono geni (righe) nella matrice x method = None
  # con un matrice di queste dimensioni infatti non è possibile applicare le procedure come TMM
  if (nrow(x) == 0 || nsamples == 1) 
    method = "none"
  
  
  # Creare uno switch, ovvero a seconda del metodo specificato intraprendere l'azione
  # Switch eliminato in questa funzione perchè voglio usare solo il TMM
  
  

  
  # Inizializzare un vettore f, tanti NA quanti il numero di samples
  f <- rep_len(NA_real_, nsamples)
  
  # Calcolare i valore di f per ogni sample.
  # Viene calcolato richiamando la funzione calcFactorTMM
  for (i in 1:nsamples)
    f[i] <- calcFactorTMM.NEW(obs = x[,i], # Selezionare l'iesimo sample/colonna
                              ref = refColumn, # Passare la colonna/sample di riferimento
                              libsize.obs = lib.size[i],  # Passare la lib size dell'iesimo sample 
                              libsize.ref = lib.size.refColumn, # Passare la libsize del sample di riferimento
                              logratioTrim = logratioTrim,  # Parametro
                              sumTrim = sumTrim, # Parametro
                              doWeighting = doWeighting,# Parametro
                              Acutoff = Acutoff # Parametro
    )
  
  ff <- f 
  f <- f/divisore.normalizzazione
  names(f) <- colnames(x)
  
  
  

  
  output <- list(
    counts = object,
    
    fattori.normalizzazione = f,
    fattore.originale = ff,
    divisore = divisore.normalizzazione,
    
    sample.riferimento = refColumn,
    libsize.sample.riferimento = lib.size.refColumn,
    
    lib.size = lib.size,
    
    dim = dim(x),
    filtered.counts = x,
    filtro = allzero
  )
  
  return(output)
}






normalizzare <- function(counts,parametri.TMM){
  # Geni comuni ad i due dataset, escludere quelli non comunis
  common.genes <- intersect(rownames(counts),
                            rownames(parametri.TMM$sample.riferimento)
  )
  counts <- counts[common.genes,]
  
  
  nuovi.parametri <- TMMbyHand.forAnotherDataset(counts[,1,drop=FALSE],
                                                 refColumn=parametri.TMM$sample.riferimento,
                                                 divisore.normalizzazione=parametri.TMM$divisore)
  
  
  # Inizializzare il dataset con la prima iterazione, in questo modo sarà più facile aggiungere i risultati
  oggettoDGEList <- DGEList(counts = counts[,1,drop = FALSE],
                            lib.size = nuovi.parametri$lib.size,
                            norm.factors = nuovi.parametri$fattori.normalizzazione
  )
  
  
  normalized_counts.test <- cpm(oggettoDGEList,normalized.lib.sizes = T,prior.count= 1 ,log = T)
  normalized_counts.test <- as.data.frame(normalized_counts.test)
  
  for (i in 2:ncol(counts)){
    # Estrarre la iesima riga di TEST e farne un oggetto edgeR
    parametri <- TMMbyHand.forAnotherDataset(counts[,i,drop=FALSE],
                                             refColumn=parametri.TMM$sample.riferimento,
                                             divisore.normalizzazione=parametri.TMM$divisore)
    
    
    oggettoDGEList <- DGEList(counts = counts[,i,drop = FALSE],
                              lib.size = nuovi.parametri$lib.size,
                              norm.factors = nuovi.parametri$fattori.normalizzazione
    )
    
    ## drop = FALSE permette di continuare a considerare l'oggetto come un dataframe e quindi mantenere i nomi di riga
    
    # Normalizzare la iesima riga di test
    normalized_counts.test_row <- cpm(oggettoDGEList,normalized.lib.sizes = T,prior.count= 1 ,log = T)
    
    # Rendere la iesima riga normalizzata un dataframe
    normalized_counts.test_row <- as.data.frame(normalized_counts.test_row)
    
    # Salvare la iesima riga
    normalized_counts.test <- cbind(normalized_counts.test, normalized_counts.test_row) 
    
  }
  return(normalized_counts.test)
}





