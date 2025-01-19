import pandas as pd 
from sklearn.model_selection import LeaveOneOut
from sklearn.svm import SVC

import sys 
#sys.path.append('...') # Add a path for the Perfomance_metrics file
from Perfomance_metrics import metrics_binary


# MODELLO SVC:  traning con LeaveOneOut, no hypertuning
def svc_loocv(train_X, # Variabili esplicative di TRAIN 
              train_y, # Variabili target di TRAIN 

              test_X,  # Variabili esplicative di TEST 
              test_y,  # Variabili target di TEST 
              random_state: int = 0
              ) -> dict:
   """  
   Funzione prende dati di train e test ed addestra e testa un modello SVC.  
   - 2. Addestra il modello SVC sui dati di train
   - 3. Testa il modello sui dati di test ricavando 

   Input:
      - train_X: Variabili esplicative di TRAIN 
      - train_y: Variabili target di TRAIN 

      - test_X:  Variabili esplicative di TEST 
      - test_y:  Variabili target di TEST 

      - random_state: random_state
   """
   # region Controllo input  
   if not isinstance(random_state,int):
      raise TypeError("random_state deve essere un numero decimale")

   if not isinstance(labels,(dict,type(None))):
      raise TypeError("labels deve essere un dizionario")
   # endregion  

   # region Training  
   # Definisci LOOCV
   loo = LeaveOneOut()
   loo.get_n_splits(train_X)


   # Definisci le liste di predizioni e valori reali
   y_true, y_pred = [], []

   #  Esegui il ciclo LOOCV per vedere come  il modello performa sui dati di training
   for train_index, test_index in loo.split(train_X):
      X_train, X_test = train_X[train_index], train_X[test_index]  
      y_train, y_test = train_y[train_index], train_y[test_index]  

      # Crea il modello con i migliori parametri trovati
      model = SVC( random_state=random_state)

      # Allena il modello
      model.fit(X_train, y_train)

      # Predizione
      yhat = model.predict(X_test)

      # Salva le predizioni e i valori reali
      y_true.append(y_test[0]) # Valore reale
      y_pred.append(yhat[0])   # Valore predetto


   metriche_training = metrics_binary(test_y=y_true, y_pred=y_pred)

   # endregion

   # region Test  
   # Definire il modello per l'evaluation

   final_model = SVC( random_state = random_state,probability=True)

   # Addestrare il modello finale su tutti i dati
   final_model.fit(train_X, train_y)

   y_pred = final_model.predict(test_X)

   # Evaluate the performance metrics 
   metriche_test= metrics_binary(test_y=y_true, y_pred=y_pred)
   # endregion

   output = {"performance_training":metriche_training,
            "performance_test":metriche_test,
            "final_model":final_model
         }

   # Ritorna le metriche finali
   return output
