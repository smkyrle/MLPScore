
# random forest for feature importance on a regression problem
from sklearn.datasets import make_regression
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import RFECV
import matplotlib.pyplot as plt
import pandas as pd
import sys
import numpy as np
from itertools import repeat
# define dataset
def importance(data, steps, RFjobs, RFECVjobs, cross_val, destination):

	df = pd.read_hdf(data, key='df', mode='r')
	label_headers = ['PDBCode','BindingDataType', 'BindingValue','BindingUnits','Label', 'Database', 'Binding_Data_in_uM']
	X = df.drop(label_headers , axis=1)
	headers = X.columns
	y = df['Label'].copy()
	# define the model
	model = RandomForestClassifier(n_estimators=50, n_jobs = RFjobs)
	selector = RFECV(model, step=steps, cv=cross_val, n_jobs=RFECVjobs, scoring='roc_auc', verbose=2)
	selector = selector.fit(X, y)
	print(f"Optimal number of features : {selector.n_features_}")
	choices = list(selector.support_)
	ranks = list(selector.ranking_)
	results = pd.DataFrame({'Feature':headers, 'Choices':choices, 'Rank':ranks})
	results.to_csv(f'{destination}RFECVFeatureImportances.csv')

	min_features_to_select = 1
	score = list(selector.grid_scores_)
	num_feature = np.linspace(min_features_to_select, len(headers), len(score))
	scores = pd.DataFrame({'Feature Number': num_feature, 'ROC AUC Scores':score})
	scores.to_csv(f'{destination}ROC_AUC_Scores.csv')

	graph(selector, headers, destination)

def graph(selector, headers, destination):
	min_features_to_select = 1
	# save accuracies of feature counts to dataframe
	plt.figure()
	plt.xlabel("Number of features selected")
	plt.ylabel("Cross validation score (nb of correct classifications)")
	plt.plot(np.linspace(min_features_to_select, len(headers), len(selector.grid_scores_)),
	         selector.grid_scores_)
	print(selector.grid_scores_)
	print(len(selector.grid_scores_))
	plt.savefig(f'{destination}RFECVFeatureImportances.png')

def parse_args(args):

	data = args[args.index('-loc') + 1]

	destination = args[args.index('-des') + 1]

	RFjobs = int(args[args.index('-RFjobs') + 1])

	RFECVjobs = int(args[args.index('-RFECVjobs') + 1])

	steps = int(args[args.index('-step') + 1])

	cross_val = int(args[args.index('-cv') + 1])

	return data, destination, RFjobs, RFECVjobs, steps,cross_val

def main():

	data, destination, RFjobs, RFECVjobs, steps, cross_val = parse_args(sys.argv)

	importance(data, steps, RFjobs, RFECVjobs, cross_val, destination)


if __name__ == '__main__':
	main()
