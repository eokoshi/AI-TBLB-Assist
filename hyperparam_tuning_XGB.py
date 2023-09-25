import time

start = time.time()
import xgboost as xgb
from sklearn.model_selection import GridSearchCV
import pandas as pd

data = pd.read_excel("inputdata.xlsx")
data = data[data.target != 3]
data["maligvsbenign"] = [1 if x < 3 else 0 for x in data.target]

X = data[ [ "intbronch", "plasmacellinfil", "eosinoinfil", "lymphoidagg", "fibroelastosis", "op",] ]
y = data["maligvsbenign"]
XGB_params = {
    "n_estimators": [10, 50, 100, 1000],
    "max_depth": [0, 3, 5, 10],
    "gamma": [0, 0.0001, 0.001, 0.01, 0.1],
    "reg_alpha": [1e-2, 0.1, 1, 10, 100],
    "reg_lambda": [1e-2, 0.1, 1, 10, 100],
    "learning_rate": [0.0001, 0.001, 0.01, 0.1,0.3,0.5, 1],
    "min_child_weight":[1, 5, 10, 100],
    "subsample":[0.5,0.7,1],
    "tree_method":["exact", "hist"],
    "objective":["binary:logistic","binary:logitraw","binary:hinge"],
}

classifier = xgb.XGBClassifier(n_jobs=1, seed=84, booster="gbtree")
gs = GridSearchCV(
    classifier,
    XGB_params,
    scoring=("roc_auc", "accuracy"),
    n_jobs=-1,
    cv=5,
    verbose=10,
    refit=False,
)
gs_fitted = gs.fit(X, y)
results = pd.DataFrame(gs_fitted.cv_results_)
results.to_csv(f"gsresults_xgb.csv")
print(time.time() - start)