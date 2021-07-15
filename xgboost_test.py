import numpy as np
import pandas as pd
import xgboost as xgb
from joblib import load
from xgboost.sklearn import XGBClassifier
from sklearn.metrics import precision_recall_fscore_support, matthews_corrcoef

print(xgb.__version__) # need to be '1.3.1'

with open("../final_selected_probes.txt", 'r') as f:
    selected_probes = f.readlines()
selected_probes = [probe.strip() for probe in selected_probes]

# load XGBoost model
xgb_model = load(open("../xgb_model_8.model", 'rb'))

#fn = "GSE38266" # 485577×42 -> 276129×42
#fn = "GSE77955"
#fn = "GSE52865" # 485512×57 -> 276129×57
#fn = "GSE72874-GPL13534" # 373561×250 -> 237020×200
#fn = "GSE61441" # 229844×92 -> 157528×92
#fn = "GSE75041" # 485577×66 -> 276129×66
fn = "GSE97466" # 449351×141 -> 263522×141
test_x = np.load(fn+'_m.npy')
test_x = pd.DataFrame(data=np.transpose(test_x), columns=selected_probes)

true_label = ['BLCA-Normal', 'BLCA-Tumor', 'BRCA-Normal', 'BRCA-Tumor',
              'COAD-Normal', 'COAD-Tumor', 'ESCA-Normal', 'ESCA-Tumor',
              'HNSC-Normal', 'HNSC-Tumor', 'KIRC-Normal', 'KIRC-Tumor',
              'KIRP-Normal', 'KIRP-Tumor', 'LIHC-Normal', 'LIHC-Tumor',
              'LUAD-Normal', 'LUAD-Tumor', 'LUSC-Normal', 'LUSC-Tumor',
              'PRAD-Normal', 'PRAD-Tumor', 'THCA-Normal', 'THCA-Tumor',
              'UCEC-Normal', 'UCEC-Tumor']

pred = xgb_model.predict(test_x)

with open("label_"+fn+".txt", "r") as f:
    label = f.readlines()
label = [i.strip() for i in label]

print(precision_recall_fscore_support(true, pred, average=None, labels=np.unique(pred))[:3])