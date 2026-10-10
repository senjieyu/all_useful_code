import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.svm import SVR
from sklearn.metrics import mean_squared_error
from deap import base, creator, tools, algorithms

# 适应度函数，评估模型性能
def evaluate_individual(individual, X_train, y_train, X_val, y_val):
    selected_features = [i for i, val in enumerate(individual) if val]
    if not selected_features:
        return float('inf'),
    
    clf = SVR(kernel='rbf')
    clf.fit(X_train[:, selected_features], y_train)
    y_pred = clf.predict(X_val[:, selected_features])
    mse = mean_squared_error(y_val, y_pred)
    return mse,

# 遗传算法设置
creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
creator.create("Individual", list, fitness=creator.FitnessMin)

toolbox = base.Toolbox()
toolbox.register("attr_bool", np.random.choice, [True, False])

# 从CSV文件读取数据
data = pd.read_csv('selected_15_data.csv')
y = data.iloc[:, 1].values
X = data.iloc[:, 2:].values
num_features = X.shape[1]

toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_bool, n=num_features)
toolbox.register("population", tools.initRepeat, list, toolbox.individual)

# 准备训练和验证数据
X_train, X_val, y_train, y_val = train_test_split(X, y, test_size=0.2, random_state=42)

toolbox.register("evaluate", evaluate_individual, X_train=X_train, y_train=y_train, X_val=X_val, y_val=y_val)
toolbox.register("mate", tools.cxTwoPoint)
toolbox.register("mutate", tools.mutFlipBit, indpb=0.05)
toolbox.register("select", tools.selTournament, tournsize=3)

def main():
    population_size = 50
    num_generations = 20
    
    population = toolbox.population(n=population_size)
    algorithms.eaSimple(population, toolbox, cxpb=0.5, mutpb=0.2, ngen=num_generations, verbose=True)
    
    best_individual = tools.selBest(population, k=1)[0]
    selected_features = [i for i, val in enumerate(best_individual) if val]
    
    print("GA Selected features:", selected_features)
    
    # 保存筛选后的特征到CSV文件
    selected_data = data.iloc[:, [0, 1] + selected_features]
    selected_data.to_csv('selected_features.csv', index=False)
    print("Selected features saved to 'selected_features.csv'")
    
    # 在最终选择的特征上训练SVM模型
    clf = SVR(kernel='rbf')
    clf.fit(X_train[:, selected_features], y_train)
    y_pred = clf.predict(X_val[:, selected_features])
    mse = mean_squared_error(y_val, y_pred)
    print("Validation MSE:", mse)
    
if __name__ == "__main__":
    main()
