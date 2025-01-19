import numpy as np
from sklearn.metrics import confusion_matrix, accuracy_score, f1_score, recall_score,precision_score

def metrics_binary(test_y, y_pred) -> dict: 
   """Ottenere le metriche  dalla confusion matrix preferito.
      Calcola: 
         - recall
         - specificity
         - precision
         - F1Score
         - Accuracy

      Utilizza le funzioni:
        - Non utilizza funzioni esterne, utilizza quelle dei pacchetti ad esempio sklearn (classification_report)

      Prende in input:
        - test_y: veri valori di y
        - y_pred: valori predetti di y

      Output:
         - dizionario: dizionario  con le metriche di performance
   """ 
   cm = confusion_matrix(test_y, y_pred)
   # Estrazione dei valori dalla confusion matrix
   TN, FP, FN, TP = cm.ravel()
   # Calcolare la specificità
   specificity = TN / (TN + FP)

   metrics_scikit_learn= {
        'Accuracy': accuracy_score(test_y, y_pred),
        'Precision': precision_score(test_y, y_pred),
        'Recall': recall_score(test_y, y_pred),
         'Specificity':specificity,
        'F1-Score': f1_score(test_y, y_pred)
    }

   return metrics_scikit_learn

