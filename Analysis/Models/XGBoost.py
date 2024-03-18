from configparser import ConfigParser
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import joblib
from xgboost import XGBClassifier
from sklearn.metrics import f1_score, accuracy_score
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import StratifiedKFold
from imblearn.over_sampling import KMeansSMOTE
from RandomForest import TransformTrainData
from MakeFile.FileSaver import FileSaver

def F1_eval(preds, dtrain):
    labels = dtrain.get_label()
    err = 1 - f1_score(labels, np.round(preds), average = "binary")
    return 'f1_err', err

if __name__ == "__main__":
    _configPath = "Analysis/Models/config.ini"
    _config = ConfigParser()
    _config.read(_configPath)   

    # filter out the CpG sites
    _aucDf = pd.read_csv(_config["Paths"]["AUC_GROUP_DATA_PATH"], usecols=["CpG", "DNAm"])
    _aucDf = _aucDf[_aucDf['DNAm'] == "hyper"]
    keepFeature = _aucDf["CpG"].tolist()
    keepFeature.append("cancer")


    # read the training data
    _trainDf = pd.read_csv(_config["Paths"]["TRAIN_BETA_DATA_PATH"], index_col=0)
    _trainDf = TransformTrainData(_trainDf, 50)
    _trainDf = _trainDf[_trainDf.columns.intersection(keepFeature)]
    _trainDf = _trainDf.iloc[1:]


    # read the testing data
    _testDf = pd.read_csv(_config["Paths"]["TEST_BETA_DATA_PATH"], index_col=0)
    _testDf = TransformTrainData(_testDf, 50)
    _testDf = _testDf[_testDf.columns.intersection(keepFeature)]
    _testDf = _testDf.iloc[1:]

    # split the training, testing data into X and Y
    _trainX = _trainDf.drop(columns=["cancer"])
    _trainY = _trainDf["cancer"]
    _testX = _testDf.drop(columns=["cancer"])
    _testY = _testDf["cancer"]

    # oversample the training data
    _trainX, _trainY = KMeansSMOTE().fit_resample(_trainX, _trainY)

    # train the model
    # _xgboostModel = XGBClassifier(
    #     n_estimators = 125, 
    #     max_depth = 4, 
    #     min_child_weight = 2, 
    #     subsample = 0.9,
    #     colsample_bytree = 0.8,
    #     learning_rate= 0.05,
    # )
    _xgboostModel = XGBClassifier(
        n_estimators = 125, 
        max_depth = 4, 
        min_child_weight = 2, 
        subsample = 0.9,
        colsample_bytree = 0.8,
        reg_alpha= 0.1,
        reg_lambda= 1,
        learning_rate= 0.05,
        objective= 'binary:logistic'
    )
    eval_set = [(_testX, _testY)]
    _xgboostModel.fit(_trainX, _trainY, eval_metric=F1_eval, eval_set=eval_set, verbose=True)

    trainPredicted = _xgboostModel.predict(_trainX)
    accuracy = accuracy_score(_trainY, trainPredicted)
    f1 = f1_score(_trainY, trainPredicted, average = "binary")
    print(accuracy)
    print("F1: ", f1)
    
    trainPredicted = _xgboostModel.predict(_testX)
    accuracy = accuracy_score(_testY, trainPredicted)
    f1 = f1_score(_testY, trainPredicted, average = "binary")
    print(accuracy)
    print("F1: ", f1)

    # save the model
    xgboostPath = _config.get('Paths', 'XGBOOST_PATH')
    joblib.dump(_xgboostModel, xgboostPath)

    # _importance = _xgboostModel.feature_importances_
    # _weight = _xgboostModel.get_booster().get_score(importance_type='weight')
    # _gain = _xgboostModel.get_booster().get_score(importance_type='gain')
    # _cover = _xgboostModel.get_booster().get_score(importance_type='cover')

    # df = pd.DataFrame({
    #     'CpG': list(_weight.keys()),
    #     'weight': list(_weight.values()),
    #     'gain': [_gain.get(feature, 0) for feature in _weight.keys()],
    #     'cover': [_cover.get(feature, 0) for feature in _weight.keys()]
    #     'importance': _importance
    # })


    # # Parameter Tuning
    # model = XGBClassifier(
    #     n_estimators = 125, 
    #     max_depth = 4, 
    #     min_child_weight = 2, 
    #     subsample = 0.9,
    #     colsample_bytree = 0.8,
    #     reg_alpha= 0.1,
    #     reg_lambda= 1,
    #     learning_rate= 0.05,
    #     objective= 'binary:logistic'
    # )
    # grid_search = GridSearchCV(model, scoring="f1", n_jobs=-1, cv=5)
    # grid_result = grid_search.fit(_trainX, _trainY)

    # # summarize results
    # print("results: {0}".format(grid_result.cv_results_))
    # print('best parameter: {0}'.format(grid_search.best_params_))
    # print('best score: {0}'.format(grid_search.best_score_))
