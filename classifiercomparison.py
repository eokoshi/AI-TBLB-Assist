#%% Setup

import pandas as pd
from sklearn.metrics import accuracy_score, roc_auc_score, log_loss, ConfusionMatrixDisplay, confusion_matrix
from sklearn.model_selection import train_test_split

#must save files as csv with commas as delimiters
data = pd.read_excel("inputdata.xlsx", index_col= "ID")

#splitting data into features and target
X = data.drop(["target"], axis = 1)
y = data["target"]

target_names = {
    1:"malignant",
    2:"probably malignant",
    3:"unclear",
    4:"probably benign",
    5:"benign"}

#%% Classifiers


from sklearn.neural_network import MLPClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt

names = [
    "Nearest Neighbors",
    "Linear SVM",
    "RBF SVM",
    "Gaussian Process",
    "Decision Tree",
    "Random Forest",
    "Neural Net",
    "AdaBoost",
    "Naive Bayes",
]


scores=[]

for asdf in range(2000):
	print(asdf)
	X_train, X_test, y_train, y_test = train_test_split(X,y, stratify=data["target"])
	classifiers = [
	    KNeighborsClassifier(5, algorithm="auto"),
	    SVC(kernel="linear", C=0.025, probability=True),
	    SVC(gamma=2, C=1, probability=True),
	    GaussianProcessClassifier(1.0 * RBF(1.0), n_restarts_optimizer=10),
	    DecisionTreeClassifier(max_depth=5),
	    RandomForestClassifier(max_depth=5, n_estimators=10),
	    MLPClassifier(alpha=1, max_iter=1000),
	    AdaBoostClassifier(),
	    GaussianNB(),
	]
	for name, clf in zip(names, classifiers):
		if name=="Linear SVM" or name=="Neural Net" or name=="Nearest Neighbors":
			clf = make_pipeline(StandardScaler(), clf)
		clf.fit(X_train, y_train)
		score = clf.score(X_test, y_test)
		proba = clf.predict_proba(X_test)
		auc = roc_auc_score(y_test, proba, multi_class="ovr")
		output= [name, score, auc]
		scores.append(output)

scores = pd.DataFrame(scores)

#%% Creating DF for plots

accuracy = [scores[1][i] for i in range(len(scores))]

auc = [scores[2][i] for i in range(len(scores))]
# c = []
# for i in range(len(scores)): c.append(namesplot[scores[i][0]])
namesplot = [scores[0][i] for i in range(len(scores))]

plotdf = pd.DataFrame([accuracy,auc,namesplot]).T
plotdf = plotdf.rename(columns={0:"Accuracy",1:"AUC",2:"Algorithm"})

medians = plotdf.groupby("Algorithm")[["Accuracy","AUC"]].median()
means = plotdf.groupby("Algorithm")[["Accuracy","AUC"]].mean()

#%% seaborn

import seaborn as sns
import matplotlib.pyplot as plt
sns.set_theme(style="whitegrid")


# Draw a scatter plot while assigning point colors and sizes to different
# variables in the dataset
fig, ax = plt.subplots(figsize=(8,7))
fig.set_dpi(600)
sns.despine(fig, left=True, bottom=True)
sns.scatterplot(x="Accuracy", y="AUC",
                hue="Algorithm",
                palette="tab10",
				style="Algorithm",
                linewidth=0, s=150,
                data=means, ax=ax, edgecolors="grey")