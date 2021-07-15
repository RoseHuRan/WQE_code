import numpy as np
import pandas as pd
import xgboost as xgb
from xgboost.sklearn import XGBClassifier
from sklearn import metrics
from sklearn.model_selection import GridSearchCV, cross_val_score
import random
import sys

'''
all_m = np.load('all_m.npy') #(276129, 6203)
all_m = np.transpose(all_m)

sel = VarianceThreshold(threshold=(.8 * (1 - .8)))

from numpy import inf
all_m[all_m == inf] = 10000
sel_m = sel.fit_transform(all_m)

with open('sel_m.npy', 'wb') as f:
    np.save(f, sel_m)
'''
seed = int(sys.argv[1])

with open("all_probe_names.txt", 'r') as f:
    probe_names = f.readlines()
probe_names = [i.strip() for i in probe_names]

all_m = np.load("all_m_new.npy")
all_m = np.transpose(all_m) #(6203, 276129)

all_m_var = np.var(all_m, axis=0)
sel_index = all_m_var > np.nanquantile(all_m_var, 0.5) #choose 50% most variable features
sel_m = all_m[:, sel_index]
sel_probe = np.array(probe_names)[sel_index]

encode_labels = np.load("encode_labels.npy") #(6203,)

random.seed(seed)
test_index = random.sample(range(len(encode_labels)), 1550) #1550
test_index.sort()

train_index = [i for i in range(len(encode_labels)) if i not in test_index] #4653


train_x = pd.DataFrame(data=sel_m[train_index], columns=sel_probe)
train_y = encode_labels[train_index]

test_x = pd.DataFrame(data=sel_m[test_index], columns=sel_probe)
test_y = encode_labels[test_index]

def modelMetrics(clf, train_x, train_y, isCv=True, cv_folds=5, early_stopping_rounds=5):  
    f = open("xgboost_out.log", "a")
    
    if isCv:  # cv chooses n_estimators 
        xgb_param = clf.get_xgb_params()
        xgb_param['num_class'] = 26
        xgtrain = xgb.DMatrix(train_x, label=train_y)  
        cvresult = xgb.cv(xgb_param, xgtrain, num_boost_round=clf.get_params()['n_estimators'], nfold=cv_folds,  
                          metrics='merror', early_stopping_rounds=early_stopping_rounds)
        clf.set_params(n_estimators=cvresult.shape[0])
        f.write("n_estimators: " + str(cvresult.shape[0]))

    #### Training ####  
    clf.fit(train_x, train_y, eval_metric='merror')

    #### Predicting ####
    predictions = clf.predict(test_x)  
    predprob = clf.predict_proba(test_x)

    #### Results ####  
    f.write("\nModel Report: seed "+str(seed)+"\n")
    f.write("Accuracy : %.4g\n" % metrics.accuracy_score(test_y, predictions))  
    f.write("ROC AUC Score (Train): %f\n" % metrics.roc_auc_score(test_y, predprob, multi_class='ovr'))
    f.close()
    
    predictions = [str(pred) for pred in predictions]
    with open('xgboost_predictions_'+str(seed)+'.txt', 'w') as f:
        f.write('\n'.join(predictions))
    
    #### Feature Importance ####
    feat_imp = pd.Series(clf.get_booster().get_fscore()).sort_values(ascending=False)
    feat_name = [i for i in feat_imp.index]
    
    return feat_name

clf1 = XGBClassifier(learning_rate=0.1, ##
                     n_estimators=100, ###
                     max_depth=10, #
                     min_child_weight=1, #
                     gamma=0, #
                     subsample=0.8,
                     colsample_bytree=0.8,
                     objective='multi:softmax',
                     nthread=8,
                     use_label_encoder=False)   #scale_pos_weight=99, 

feat_name = modelMetrics(clf1, train_x, train_y, isCv=False)

with open('feature_index_'+str(seed)+'.txt', 'w') as f:
    f.write('\n'.join(feat_name))