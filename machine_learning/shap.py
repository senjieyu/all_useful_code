

import pandas as pd, numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR
# from shap import KernelExplainer # Kernel explainer貌似做不出bar图，如果你看到了，你可以切换成kernel
from shap.explainers import Permutation, Exact # 我自己喜欢用的2个用来代替kernel的explainer
from shap import plots
import matplotlib.pyplot as plt
plt.rcParams["figure.dpi"] = 300

train = pd.read_excel("subtrain.xlsx")
train.iloc[:, 2:] = StandardScaler().fit_transform(train.iloc[:, 2:])
X_train = train.iloc[:, 2:]
Y_train = train.iloc[:, 1]

model_p = dict(
        C=24,
        gamma=0.6,
        epsilon=0.06,
    )

model = SVR(**model_p).fit(X_train, Y_train)

explainer = Permutation(model.predict, X_train, feature_names=X_train.columns)
shap_values = explainer(X_train)

# bar图
plots.bar(shap_values)
# 显示前n个feature
plots.bar(shap_values, max_display=12)
# 显示第一个样本的变量贡献
plots.bar(shap_values[0])

# 蜂群图
'''plots.beeswarm(shap_values)
# 显示前n个geature
plots.beeswarm(shap_values, max_display=20)
# 蜂群图可以传入其他颜色
plots.beeswarm(shap_values, color=plt.get_cmap("cool"))

# decision plot
plots.decision(shap_values.base_values[0], shap_values.values, 
               feature_names=X_train.columns.values)'''

# 哑铃图
plots.violin(shap_values.values, X_train)

# 单个样本瀑布图
'''plots.waterfall(shap_values[0,:])

# 热力图
plots.heatmap(shap_values)

# 每个样本的敏感性分析，1是feature的index
plots.partial_dependence(2, model.predict, X_train, feature_names=X_train.columns.values)
# 敏感性分析
plots.partial_dependence(1, model.predict, X_train, feature_names=X_train.columns.values, ice=False)
# index可以换成feature name
plots.partial_dependence("Cu", model.predict, X_train, feature_names=X_train.columns.values, ice=False)'''

# 每个feature的散点图，0是feature的index
plots.scatter(shap_values[:, 0])
# index可以换成feature name
plots.scatter(shap_values[:, "c6"])

# 如果你要自己作图，有几点可以注意
# shap_values分为几个部分
'''shapley_values = shap_values.values
origin_feature_values = shap_values.data
base_values = shap_values.base_values

# 如果你要作bar图，也就是需要变量重要性
importances = np.abs(shapley_values).mean(axis=0)

# 如果你需要index为0的feature的散点图，这里的index不能用feature name
plot_x = origin_feature_values[:, 0] # feature值
plot_y = shapley_values[:, 0] # feature的shapley values
fig, ax = plt.subplots(1)
ax.scatter(plot_x, plot_y)'''

# 其他图的重现比较麻烦了可能，我没研究过，你自己去研究吧