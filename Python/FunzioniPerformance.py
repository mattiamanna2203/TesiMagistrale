import pandas as pd
import numpy as np
from sklearn.metrics import accuracy_score,confusion_matrix,  accuracy_score, f1_score
def evaluation_multiclass(test_X, test_y, model):
   """Ottenere le metriche aggiuntive dalla confusion matrix preferito a  classification_report in 
      questo modo si può calcolare anche la sensitivity.
      Calcola: 
         - sensitivity (recall)
         - specificity
         - precision
         - F1Score
         - Accuracy

      Utilizza le funzioni:
        - metrics_from_confusion_matrix
        - get_overall_metrics

      Prende in input:
        - ConfusioMatrix

      Output:
         - metriche: Pandas dataframe con le metriche di performance per ogni classe
         - overall_metrics: dataframe con le metriche overall
   """ 
   # Utilizzare il modello e ottenere le predizioni su dati mai visti dal modello
   y_pred = model.predict(test_X)

   # Classi
   labels = np.unique(test_y)  

   # Confusion matrix 
   cm = confusion_matrix(test_y, y_pred, labels=labels)

   # Ottenere le metriche dalla confusion matrix
   metriche  = metrics_from_confusion_matrix(cm)
   overall_metrics = get_overall_metrics(test_y, y_pred)


   return metriche, cm, overall_metrics

def get_overall_metrics(y_true_fold, y_pred_fold):
   """ Prende in input i veri valori di y e quelli predetti e calcola la metriche generali come:
      - Accuracy overall
      - macroF1score
      - microF1score
      - weightedF1score

      Utilizza le funzioni:
        - Non utilizza funzioni esterne, utilizza quelle dei pacchetti ad esempio sklearn


      Prende in input:
         - y_true_fold: veri valori di y
         - y_pred_fold: valori predetti di y

      Output:
         - overall_metrics_dataframe: pandas dataframe contenente tutte le metriche
   """
   # Overall metrics tramite le built in functions, si spiega da solo
   overall_accuracy = accuracy_score(y_true_fold, y_pred_fold)
   macro_f1 = f1_score(y_true_fold, y_pred_fold, average="macro")
   micro_f1 = f1_score(y_true_fold, y_pred_fold, average="micro")
   weighted_f1 = f1_score(y_true_fold, y_pred_fold, average="weighted")


   # Salvare le metriche overall in un dictionary
   overall_metrics = {
      "accuracy": overall_accuracy,
      "macroF1score": macro_f1,
      "microF1score": micro_f1,
      "weightedF1score": weighted_f1,
   }

   overall_metrics_dataframe = pd.DataFrame.from_dict(overall_metrics,orient="index").T

   return overall_metrics_dataframe

def metrics_from_confusion_matrix(ConfusioMatrix):
   """Ottenere le metriche aggiuntive dalla confusion matrix preferito a  classification_report in 
      questo modo si può calcolare anche la sensitivity.
      Calcola: 
         - sensitivity (recall)
         - specificity
         - precision
         - F1Score
         - Accuracy

      Utilizza le funzioni:
        - Non utilizza funzioni esterne, utilizza quelle dei pacchetti ad esempio sklearn

      Prende in input:
        - ConfusioMatrix

      Output:
         - metriche: Pandas dataframe con le metriche di performance per ogni classe
   """ 
   # Inizializzare il dataframe ove verranno salvate 
   metriche = pd.DataFrame()

   # Total number of samples
   total_samples = np.sum(ConfusioMatrix)  

   # region Calcolare le metriche
   # Iterare per  ogni riga della confusion matrix, quindi per ogni classe 
   for i in range(len(ConfusioMatrix)):
      TP = ConfusioMatrix[i, i]  # True Positives for class i
      FP = ConfusioMatrix[:, i].sum() - TP  # False Positives for class i (other classes predicted as class i)
      FN = ConfusioMatrix[i, :].sum() - TP  # False Negatives for class i (class i predicted as other classes)
      TN = total_samples - (TP + FP + FN)  # True Negatives for class i

      # Calcolare le metriche, i nomi si spiegano da soli, inoltre mettere dei check aggiuntivi sulla divisione per zero
      specificity = TN / (TN + FP) if (TN + FP) > 0 else 0
      precision = TP / (TP + FP) if (TP + FP) > 0 else 0
      recall = TP  / (TP+FN) if (TP + FN) > 0 else 0
      f1_score = 2 * ((precision *recall)/(precision  + recall))  if (precision  + recall) > 0 else 0


      # Un check per vedere se i risultati corrispondono
      if TP+TN+FP+FN != total_samples:
         raise ValueError("Error of some kind occurred")

      # Calcolare l'accuracy
      accuracy = (TP+TN)/(total_samples)

      # Salvare le metriche ottenute per la classe i-esima
      row = pd.DataFrame({ "recall":[recall],
                           "specificity":[specificity],
                           "precision":[precision],
                           "F1Score":[f1_score],
                           "accuracy":[accuracy]})


      # Salvare  le metriche nel dataframe principale
      metriche=pd.concat([metriche,row])
   # endregion
   # Assegnare gli index al dataframe, senza questa riga erano tutti zeri.
   metriche.index = [float(i) for i in range(len(ConfusioMatrix))]

   return metriche   