import pandas as pd
import autosklearn.classification
from sklearn.model_selection import train_test_split
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import classification_report
from sklearn.svm import SVC

import sklearn.metrics

file = pd.read_excel('/Users/zhenjia/Desktop/20160818_NAI.xlsx', sheet_name='Sheet1')
x1 = file['bloom_NAI_SS490:float']
x2 = file['bbpindex555']
y = file['type']
X = pd.concat([x1, x2], axis=1)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3,
                                                    random_state=7)  # test size代表了测试集在整个数据集中的占比，这里选取的是0.3即30%


tuned_parameters = [{'kernel': ['rbf'], 'gamma': [1e-3, 1e-4, 1e-5, 1e-2, 1e-1,0.05], 'C': [1e+3,1e+4,1e+5,1e+6,1e+7]},
                    {'kernel': ['linear'], 'C': [1, 10, 100, 1000]},
                    {'kernel': ['sigmoid'], 'gamma': [1e-3, 1e-4, 1e-5, 1e-2, 1e-1]},
                    {'kernel': ['poly'], 'degree': [2, 3, 4, 5], 'gamma': [1e-3, 1e-4, 1e-5, 1e-2, 1e-1]}]
scores=['accuracy']
for score in scores:
    print("# Tuning hyper-parameters for %s" % score)
    print()

     # 调用 GridSearchCV，将 SVC(), tuned_parameters, cv=5, 还有 scoring 传递进去，
    clf = GridSearchCV(SVC(), tuned_parameters, cv=5,
                       scoring=score)
    # 用训练集训练这个学习器 clf
    clf.fit(X_train, y_train)

    print("Best parameters set found on development set:")
    print()

    # 再调用 clf.best_params_ 就能直接得到最好的参数搭配结果
    print(clf.best_params_)

    print()
    print("Grid scores on development set:")
    print()
    means = clf.cv_results_['mean_test_score']
    stds = clf.cv_results_['std_test_score']

    # 看一下具体的参数间不同数值的组合后得到的分数是多少
    for mean, std, params in zip(means, stds, clf.cv_results_['params']):
        print("%0.3f (+/-%0.03f) for %r"
              % (mean, std * 2, params))

    print()

    print("Detailed classification report:")
    print()
    print("The model is trained on the full development set.")
    print("The scores are computed on the full evaluation set.")
    print()
    y_true, y_pred = y_test, clf.predict(X_test)

    # 打印在测试集上的预测结果与真实值的分数
    print(classification_report(y_true, y_pred))

    print()
# ## Tuning hyper-parameters for accuracy
# Best parameters set found on development set:
# {'C': 10000000.0, 'gamma': 0.1, 'kernel': 'rbf'}
# Grid scores on development set:
# 0.559 (+/-0.010) for {'C': 1000.0, 'gamma': 0.001, 'kernel': 'rbf'}
# 0.559 (+/-0.010) for {'C': 1000.0, 'gamma': 0.0001, 'kernel': 'rbf'}
# 0.559 (+/-0.010) for {'C': 1000.0, 'gamma': 1e-05, 'kernel': 'rbf'}
# 0.559 (+/-0.010) for {'C': 1000.0, 'gamma': 0.01, 'kernel': 'rbf'}
# 0.619 (+/-0.125) for {'C': 1000.0, 'gamma': 0.1, 'kernel': 'rbf'}
# 0.557 (+/-0.016) for {'C': 1000.0, 'gamma': 0.05, 'kernel': 'rbf'}
# 0.559 (+/-0.010) for {'C': 10000.0, 'gamma': 0.001, 'kernel': 'rbf'}
# 0.559 (+/-0.010) for {'C': 10000.0, 'gamma': 0.0001, 'kernel': 'rbf'}
# 0.559 (+/-0.010) for {'C': 10000.0, 'gamma': 1e-05, 'kernel': 'rbf'}
# 0.559 (+/-0.010) for {'C': 10000.0, 'gamma': 0.01, 'kernel': 'rbf'}
# 0.627 (+/-0.096) for {'C': 10000.0, 'gamma': 0.1, 'kernel': 'rbf'}
# 0.606 (+/-0.101) for {'C': 10000.0, 'gamma': 0.05, 'kernel': 'rbf'}
# 0.559 (+/-0.010) for {'C': 100000.0, 'gamma': 0.001, 'kernel': 'rbf'}
# 0.559 (+/-0.010) for {'C': 100000.0, 'gamma': 0.0001, 'kernel': 'rbf'}
# 0.559 (+/-0.010) for {'C': 100000.0, 'gamma': 1e-05, 'kernel': 'rbf'}
# 0.621 (+/-0.132) for {'C': 100000.0, 'gamma': 0.01, 'kernel': 'rbf'}
# 0.672 (+/-0.089) for {'C': 100000.0, 'gamma': 0.1, 'kernel': 'rbf'}
# 0.650 (+/-0.104) for {'C': 100000.0, 'gamma': 0.05, 'kernel': 'rbf'}
# 0.584 (+/-0.057) for {'C': 1000000.0, 'gamma': 0.001, 'kernel': 'rbf'}
# 0.557 (+/-0.016) for {'C': 1000000.0, 'gamma': 0.0001, 'kernel': 'rbf'}
# 0.574 (+/-0.077) for {'C': 1000000.0, 'gamma': 1e-05, 'kernel': 'rbf'}
# 0.682 (+/-0.062) for {'C': 1000000.0, 'gamma': 0.01, 'kernel': 'rbf'}
# 0.783 (+/-0.098) for {'C': 1000000.0, 'gamma': 0.1, 'kernel': 'rbf'}
# 0.744 (+/-0.070) for {'C': 1000000.0, 'gamma': 0.05, 'kernel': 'rbf'}
# 0.629 (+/-0.098) for {'C': 10000000.0, 'gamma': 0.001, 'kernel': 'rbf'}
# 0.588 (+/-0.087) for {'C': 10000000.0, 'gamma': 0.0001, 'kernel': 'rbf'}
# 0.591 (+/-0.180) for {'C': 10000000.0, 'gamma': 1e-05, 'kernel': 'rbf'}
# 0.766 (+/-0.096) for {'C': 10000000.0, 'gamma': 0.01, 'kernel': 'rbf'}
# 0.815 (+/-0.077) for {'C': 10000000.0, 'gamma': 0.1, 'kernel': 'rbf'}
# 0.791 (+/-0.095) for {'C': 10000000.0, 'gamma': 0.05, 'kernel': 'rbf'}
# 0.559 (+/-0.010) for {'C': 1, 'kernel': 'linear'}
# 0.559 (+/-0.010) for {'C': 10, 'kernel': 'linear'}
# 0.559 (+/-0.010) for {'C': 100, 'kernel': 'linear'}
# 0.559 (+/-0.010) for {'C': 1000, 'kernel': 'linear'}
# 0.559 (+/-0.010) for {'gamma': 0.001, 'kernel': 'sigmoid'}
# 0.559 (+/-0.010) for {'gamma': 0.0001, 'kernel': 'sigmoid'}
# 0.559 (+/-0.010) for {'gamma': 1e-05, 'kernel': 'sigmoid'}
# 0.559 (+/-0.010) for {'gamma': 0.01, 'kernel': 'sigmoid'}
# 0.559 (+/-0.010) for {'gamma': 0.1, 'kernel': 'sigmoid'}
# 0.559 (+/-0.010) for {'degree': 2, 'gamma': 0.001, 'kernel': 'poly'}
# 0.559 (+/-0.010) for {'degree': 2, 'gamma': 0.0001, 'kernel': 'poly'}
# 0.559 (+/-0.010) for {'degree': 2, 'gamma': 1e-05, 'kernel': 'poly'}
# 0.559 (+/-0.010) for {'degree': 2, 'gamma': 0.01, 'kernel': 'poly'}
# 0.559 (+/-0.010) for {'degree': 2, 'gamma': 0.1, 'kernel': 'poly'}
# 0.559 (+/-0.010) for {'degree': 3, 'gamma': 0.001, 'kernel': 'poly'}
# 0.559 (+/-0.010) for {'degree': 3, 'gamma': 0.0001, 'kernel': 'poly'}
# 0.559 (+/-0.010) for {'degree': 3, 'gamma': 1e-05, 'kernel': 'poly'}
# 0.559 (+/-0.010) for {'degree': 3, 'gamma': 0.01, 'kernel': 'poly'}
# 0.559 (+/-0.010) for {'degree': 3, 'gamma': 0.1, 'kernel': 'poly'}
# 0.559 (+/-0.010) for {'degree': 4, 'gamma': 0.001, 'kernel': 'poly'}
# 0.559 (+/-0.010) for {'degree': 4, 'gamma': 0.0001, 'kernel': 'poly'}
# 0.559 (+/-0.010) for {'degree': 4, 'gamma': 1e-05, 'kernel': 'poly'}
# 0.559 (+/-0.010) for {'degree': 4, 'gamma': 0.01, 'kernel': 'poly'}
# 0.559 (+/-0.010) for {'degree': 4, 'gamma': 0.1, 'kernel': 'poly'}
# 0.559 (+/-0.010) for {'degree': 5, 'gamma': 0.001, 'kernel': 'poly'}
# 0.559 (+/-0.010) for {'degree': 5, 'gamma': 0.0001, 'kernel': 'poly'}
# 0.559 (+/-0.010) for {'degree': 5, 'gamma': 1e-05, 'kernel': 'poly'}
# 0.559 (+/-0.010) for {'degree': 5, 'gamma': 0.01, 'kernel': 'poly'}
# 0.559 (+/-0.010) for {'degree': 5, 'gamma': 0.1, 'kernel': 'poly'}
# Detailed classification report:
# The model is trained on the full development set.
# The scores are computed on the full evaluation set.
#               precision    recall  f1-score   support
#            0       0.75      0.81      0.78       102
#            1       0.79      0.73      0.76        99
#     accuracy                           0.77       201
#    macro avg       0.77      0.77      0.77       201
# weighted avg       0.77      0.77      0.77       201

