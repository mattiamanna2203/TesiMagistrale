import pandas as pd # Lavorare con dataframe
from sklearn.model_selection import KFold, LeaveOneOut, GridSearchCV, StratifiedKFold, train_test_split
from sklearn.metrics import accuracy_score,confusion_matrix,make_scorer, recall_score
from sklearn.svm import SVC

# Pacchetto per poter aggiungere la repository
import sys 
# Aggiungere un path file ai pacchetti così da importare le funzioni da me create
sys.path.append('/Users/mattia/Desktop/Università/Data Science in Python/__pacchetti__')

from get_performance_metrics import evaluation_multiclass,metrics_from_confusion_matrix,get_overall_metrics
# Funzioni di supporto 
# Funzione per assegnare le labels alle classi
def __set_labels__(dataframe:pd.DataFrame(),
               labels : dict
               ) -> pd.DataFrame():
   """Funzione per assegnare le labels alle classi 
      Prende un dataframe che illustra le metriche di performance di ogni classe di una classificazione
      Restituisce il dataframe con i nomi delle classi al posto delle metriche.
   """
   df = dataframe.copy() # Per evitare errori di sovrascrittura

   for index,patologia in labels.items():
      if df.loc[index].name == float(index):
         df.at[index,"index"] = patologia

   

   df.set_index("index",inplace=True,drop=True)
   df.index.name = None

   return df 

#------------#------------#------------#------------#------------#------------#------------#------------#------------#
# MODELLO SVC: Fatto HyperTuning con KFold, Testato su dati di traning con LeaveOneOut
def svc_loocv_hypKF(train_X, # Variabili esplicative di TRAIN 
                    train_y, # Variabili target di TRAIN 

                    test_X,  # Variabili esplicative di TEST 
                    test_y,  # Variabili target di TEST 

                    labels : dict = None,
                    ottimizza: str  = 'accuracy',
                    random_state: int = 0,
                    param_grid: dict  = { 'C': [0.1, 1, 10, 100],       # Esempio di valori per il parametro C
                                            'kernel': ['linear', 'rbf'],   # Sperimenta sia kernel lineari che RBF
                                            'gamma': ['scale', 'auto']    # Opzioni per il parametro gamma
                                        }
                        ):
    """  
    Funzione prende dati di train e test ed addestra e testa un modello SVC.  
    - 1. Performa KFold per tuning dei parametri su dati di train (Metodo GridSearchCV)
    - 2. Addestra il modello SVC sui dati di train
    - 3. Testa il modello sui dati di test ricavando 

    Input:
        - train_X: Variabili esplicative di TRAIN 
        - train_y: Variabili target di TRAIN 

        - test_X:  Variabili esplicative di TEST 
        - test_y:  Variabili target di TEST 

        - labels: dizionario contenente i nomi delle classi che si stanno classificando 
        - ottimizza (default accuracy) : Metrica di performance da ottimizzare durante il tuning degli hyperparameters
        - random_state: random_state
        - param_grid: Opzioni dei parametri sui quali applicare la grid search (tuning hyperparameters)  
    """
    # region Controllo input  
    if not isinstance(ottimizza,str):
        raise TypeError("tuning_score deve essere una stringa")

    if not isinstance(random_state,int):
        raise TypeError("random_state deve essere un numero decimale")

    if not isinstance(param_grid,dict):
        raise TypeError("param_grid deve essere un dizionario")

    if not isinstance(labels,(dict,type(None))):
        raise TypeError("labels deve essere un dizionario")
    # endregion  

    # region Tuning HyperParameters    
    # Definizione del KFold con 10 fold



    # Definire il Kfold da utilizzare per il grid search
    kf = StratifiedKFold(n_splits=10, shuffle=True, random_state = random_state)




    # Definizione del modello di base
    base_model = SVC(random_state = random_state)



    # Dividere i dati  per training e hypertunig 
    train_X, X_val, train_y, y_val = train_test_split(train_X,
                                                      train_y,
                                                      test_size=0.30,
                                                      random_state=random_state,
                                                      stratify=train_y)


    nsamples_train = pd.DataFrame(pd.DataFrame(train_y).value_counts())
    nsamples_train.reset_index(inplace=True,drop=True)
    nsamples_train.index = [labels[int(index)] for index in list(nsamples_train.index)]
    nsamples_train = nsamples_train.rename(columns={"count":"# samples train"})

    nsamples_val = pd.DataFrame(pd.DataFrame(y_val).value_counts())
    nsamples_val.reset_index(inplace=True,drop=True)
    nsamples_val.index = [labels[int(index)] for index in list(nsamples_val.index)]
    nsamples_val = nsamples_val.rename(columns={"count":"# samples hypertuning"})


    risultati_tr = pd.concat([nsamples_train,nsamples_val],axis=1)
  
        



    # Ottimizzare il recall della classe cancer
    custom_scorer = make_scorer(recall_score, average=None, labels=[1])


    


    # GridSearchCV per il tuning dei parametri
    #grid_search = GridSearchCV(estimator=base_model, param_grid=param_grid, cv=kf, n_jobs=-1, scoring = custom_scorer)
    grid_search = GridSearchCV(estimator=base_model, param_grid=param_grid, cv=kf, n_jobs=-1, scoring = ottimizza)
    

    
    # Esegui la ricerca sulla griglia di parametri con i dati di train
    grid_search.fit(X_val, y_val)
    # endregion  

    # region Training  
    # Definisci LOOCV
    loo = LeaveOneOut()
    loo.get_n_splits(train_X)
   

    # Definisci le liste di predizioni e valori reali
    y_true, y_pred = [], []
  
  # Esegui il ciclo LOOCV per vedere come  il modello performa sui dati di training
    for train_index, test_index in loo.split(train_X):
        X_train, X_test = train_X[train_index], train_X[test_index]  
        y_train, y_test = train_y[train_index], train_y[test_index]  

        # Crea il modello con i migliori parametri trovati
        model = SVC(**grid_search.best_params_, random_state=random_state)

        # Allena il modello
        model.fit(X_train, y_train)

        # Predizione
        yhat = model.predict(X_test)

        # Salva le predizioni e i valori reali
        y_true.append(y_test[0]) # Valore reale
        y_pred.append(yhat[0])   # Valore predetto
  
    # Calcola le metriche per il training
    confusion_matrix_training = confusion_matrix(y_true, y_pred)
    metriche_training  = metrics_from_confusion_matrix(confusion_matrix_training)

    overall_metrics_trainig = get_overall_metrics(y_true, y_pred)

    # endregion

    # region Test  
    # Definire il modello per l'evaluation
    final_model = SVC(**grid_search.best_params_, random_state = random_state)
    #final_model = SVC( random_state = random_state)
    # Addestrare il modello finale su tutti i dati
    final_model.fit(train_X, train_y)

    # Calcolare le metriche per l'evaluation
    metriche_test, confusion_matrix_test, overall_metrics_test = evaluation_multiclass(test_X, test_y, final_model) 
    
    # endregion

    # Labels, se sono state specificate utilizzarle
    if isinstance(labels,dict):
        metriche_training = __set_labels__(metriche_training,labels)
        metriche_test = __set_labels__(metriche_test,labels)


    parametri = pd.DataFrame.from_dict(data=grid_search.best_params_,orient="index").T

    output = {"performance_training":metriche_training,
              "performance_test":metriche_test,

              "samples":   risultati_tr,
    
              "overall_metrics_training":overall_metrics_trainig,
              "overall_metrics_test":overall_metrics_test,  

              "parametri":  parametri,
              "confusion_matrix_training":confusion_matrix_training,
              "confusion_matrix_test":confusion_matrix_test
            }

    # Ritorna le metriche finali
    return output

