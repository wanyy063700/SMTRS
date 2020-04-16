#!/usr/bin/env python

from __future__ import absolute_import, division, print_function

import pandas as pd
import os 
import re 
import glob
import argparse as arg
import sys
import source
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import sklearn


from pandas.plotting import scatter_matrix 
import matplotlib.pyplot as plt
from sklearn import model_selection
from sklearn.metrics import classification_report
from sklearn.metrics import confusion_matrix
from sklearn.metrics import accuracy_score
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.naive_bayes import GaussianNB
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier, forest

from sklearn.model_selection import GridSearchCV

import pickle


class fingerprint():
    def __init__(self, fp_fun, name):
        self.fp_fun = fp_fun
        self.name = name
        self.x = []

class MODEL_SELECTION():
    def __init__(self):
        self.models= []
        self.models.append(('LR', LogisticRegression()))
        self.models.append(('LDA', LinearDiscriminantAnalysis()))
        self.models.append(('KNN', KNeighborsClassifier()))
        self.models.append(('CART', DecisionTreeClassifier()))
        self.models.append(('NB', GaussianNB()))
        self.models.append(('SVM', SVC()))
        self.models.append(('RF', RandomForestClassifier()))

    def load(self, data='', labelcol=1, othercol=[0]):
        #labelcol  col index that store the label
        #othercol: col indexes that store data other than feature data and label data
        try:
            t= pd.read_csv(data)
        except:
            t= pd.read_csv(data, encoding='ISO-8859-1')
        
        columns = list(t.columns)
        label_col_name = columns[labelcol] 
        other_col_names = [ columns[i] for i in othercol]
        val_cols = []
        val_cols = [ i for i in columns if i not in ([label_col_name]+other_col_names)]
        #print (label_col_name, other_col_names)
        #print (len(val_cols), "\n", val_cols)
        self.x = t[val_cols].values
        self.y = t[label_col_name].values
            
    def cv(self, n_splits=10, seed=1000, scoring='accuracy'):
        #scoring:  https://scikit-learn.org/stable/modules/model_evaluation.html#scoring-parameter
        X_train=self.x
        Y_train = self.y
        results = []
        names = []
        results_summary=[]
        for name, model in self.models:
            kfold = model_selection.KFold(n_splits=n_splits, random_state=seed)
            cv_results = model_selection.cross_val_score(model, X_train, Y_train, cv=kfold, scoring=scoring)
            results.append(cv_results)
            names.append(name)
            msg = "%s: %f (%f)" % (name, cv_results.mean(), cv_results.std())
            print(msg)
            results_summary.append([name, scoring, cv_results.mean(), cv_results.std()])

        self.t_rst = pd.DataFrame(results_summary, columns=['Model','Scoring', 'Mean', 'STD'] )
        self.results = results
        self.model_names = names

    def export(self, filename='model_selection_result.csv'):
        self.t_rst.to_csv(filename,index=False)

class GRID_SEARCH():
    def __init__(self, model=RandomForestClassifier(), scoring='recall'):
        self.model = model
        self.score = scoring
        self.param_grid = { 
            'n_estimators': [200, 300, 400, 500],
            'max_features': [None, 'auto', 'sqrt', 'log2'],
            'max_depth' : [None, 7,8,9],
            'min_samples_split' : [2,4,6],
            'criterion' :['gini', 'entropy']
        }

        '''
        self.param_grid = { 
            'n_estimators': [200],
            'max_features': [None, 'auto'],
            'max_depth' : [None,9],
            'criterion' :['gini', 'entropy']
        }
        '''

    def load(self, data='', labelcol=1, othercol=[0]):
        #labelcol  col index that store the label
        #othercol: col indexes that store data other than feature data and label data
        try:
            t= pd.read_csv(data)
        except:
            t= pd.read_csv(data, encoding='ISO-8859-1')
        
        columns = list(t.columns)
        label_col_name = columns[labelcol] 
        other_col_names = [ columns[i] for i in othercol]
        val_cols = []
        val_cols = [ i for i in columns if i not in ([label_col_name]+other_col_names)]
        #print (label_col_name, other_col_names)
        #print (len(val_cols), "\n", val_cols)
        self.x = t[val_cols].values
        self.y = t[label_col_name].values

    def search(self):
        #scores = ['recall', 'precision']
        scores = ['recall']

        result_list = []
        for score in scores: 
            clf = GridSearchCV(estimator=self.model, param_grid=self.param_grid, cv=5, scoring='%s_macro' % score)
            clf.fit(self.x, self.y)
            print("Best parameters set found on development set:")
            print(clf.best_params_)
            print("Grid scores on development set:")
            means = clf.cv_results_['mean_test_score']
            stds = clf.cv_results_['std_test_score']
            for mean, std, params in zip(means, stds, clf.cv_results_['params']):
                print("%s: %0.3f (+/-%0.03f) for %r" % (score, mean, std * 2, params))
                result_list.append([score, mean, std, "%r" % params])
        self.t_rst = pd.DataFrame(result_list, columns=['Score', 'Mean', 'STD', 'params'] )

    def export(self, filename='grid_search_result.csv'):
        self.t_rst.to_csv(filename,index=False)


class RFML():
    def __init__(self):
        #{'criterion': 'gini', 'max_depth': None, 'max_features': 'auto', 'min_samples_split': 4, 'n_estimators': 500}
        print ("place holders")

    def load(self, data='', labelcol=1, othercol=[0]):
        #labelcol  col index that store the label
        #othercol: col indexes that store data other than feature data and label data
        try:
            t= pd.read_csv(data)
        except:
            t= pd.read_csv(data, encoding='ISO-8859-1')
        
        columns = list(t.columns)
        label_col_name = columns[labelcol] 
        other_col_names = [ columns[i] for i in othercol]
        val_cols = []
        val_cols = [ i for i in columns if i not in ([label_col_name]+other_col_names)]
        #print (label_col_name, other_col_names)
        #print (len(val_cols), "\n", val_cols)
        self.x = t[val_cols].values
        self.y = t[label_col_name].values
        self.other_col = t[other_col_names[0]].values

    def train(self):
        #rdkit
        #self.param = { 
        #    'n_estimators': 500,
        #    'max_features': 'auto',
        #    'max_depth' : None,
        #    'min_samples_split' : 2,
        #    'criterion' :'gini'
        #}
        #morgan
        self.param = { 
            'n_estimators': 200,
            'max_features': 'auto',
            'max_depth' : None,
            'min_samples_split' : 2,
            'criterion' :'entropy'
        }
        self.rf = RandomForestClassifier(n_estimators=self.param['n_estimators'], max_depth=self.param['max_depth'], min_samples_split=self.param['min_samples_split'], max_features=self.param['max_features'], criterion=self.param['criterion'], random_state=100)
        self.rf.fit(self.x, self.y)

    def load_model(self, model=''):
        if model == '':
            print ("No model file provided, program exit..."); exit()
        with open(model, 'rb') as f:
            self.rf = pickle.load(f)
        print (self.rf)

    def pred(self, testset=True):
        y_pred = self.rf.predict(self.x)
        y_pred_proba = self.rf.predict_proba(self.x)
        if testset==True:
            score = self.rf.score(self.x, self.y)
            print (score)
        p_vs_t = list(zip(self.other_col, self.y, y_pred, y_pred_proba[:,1]))
        self.t_pvt = pd.DataFrame(p_vs_t, columns=['Pair', 'True', 'Pred', 'Proba'] )


        self.t_pvt.sort_values(['Proba'], ascending=[0], inplace=True)
        print (self.t_pvt)
        if testset==True:
            print (classification_report(self.y, y_pred)) 
            print (classification_report(self.t_pvt[0:100]['True'].values,self.t_pvt[0:100]['Pred'].values )) 
            print (classification_report(self.t_pvt[0:50]['True'].values,self.t_pvt[0:50]['Pred'].values )) 
            print (classification_report(self.t_pvt[0:20]['True'].values,self.t_pvt[0:20]['Pred'].values )) 

    def export_result(self, filename='grid_search_result.csv'):
        self.t_pvt.to_csv(filename,index=False)

    def export(self, filename='model.pkl'):
        with open(filename, 'wb') as f:
            pickle.dump(self.rf, f, pickle.HIGHEST_PROTOCOL)

if __name__ == '__main__':
    opt=arg.ArgumentParser(description='Create dose response upload file from csv')
    opt.add_argument('-i','--input_file', help='input file')
    opt.add_argument('--rf',help='rna feature file')
    opt.add_argument('--lf',help='ligand feature file')
    opt.add_argument('-o','--output_file', type=str, help='output file')
    opt.add_argument('--model_selection', action='store_true', help='model selection')
    opt.add_argument('--grid_search', action='store_true', help='model selection')
    opt.add_argument('--train', action='store_true', help='model training')
    opt.add_argument('--pred', action='store_true', help='model pred')
    opt.add_argument('--real_set', action='store_true', help='is the prediction made on real set? if no, program will print scores based on real label and predicted label, if yes, just prediction')
    opt.add_argument('--model', type=str, help='model file')
    opt.add_argument('-s','--scoring', default='accuracy', type=str, choices=['accuracy', 'recall', 'roc_auc', 'precision'], help='scoring function in model selection or grid search, precision for grid search, accuracy for model selection')
    args=opt.parse_args()

    if args.model_selection:
        ms = MODEL_SELECTION()
        ms.load(data=args.input_file)

        print (args.scoring, "  ......................")
        ms.cv( scoring=args.scoring)
        ms.export(filename=args.output_file)

    elif args.grid_search:
        gs = GRID_SEARCH(model=RandomForestClassifier(random_state=100), scoring=args.scoring)
        gs.load(data=args.input_file)
        gs.search()
        gs.export(filename=args.output_file)

    elif args.train:
        m = RFML()
        m.load(data=args.input_file)
        m.train()
        m.export(filename=args.output_file)

    elif args.pred:
        m = RFML()
        m.load(data=args.input_file)
        m.load_model(model=args.model)
        if args.real_set:
            m.pred(testset=False)
        else:
            m.pred()
        m.export_result(filename=args.output_file)
