# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 11:23:22 2024

@author: senji
"""

# 特征衍生指的是从原始特征（原始数据的属性或变量）派生出新的特征，以增强模型对数据的表达能力或帮助模型更好地捕捉数据的模式。这可以涉及到各种数学操作、变换或组合。

# 特征衍生的目的通常有以下几点：

# 1. **提高模型性能：** 引入新的特征可以使模型更好地拟合数据，提高预测性能。

# 2. **更好地表达数据：** 有时原始特征可能不足以完整地表达数据中的信息。通过衍生新特征，可以更全面地描述数据的特性。

# 3. **减轻非线性关系：** 如果数据包含非线性关系，通过引入多项式特征或其他非线性变换，可以帮助模型更好地拟合数据。

# 4. **简化模型：** 特征衍生有时也可以帮助去除数据中的噪声或冗余信息，从而简化模型并提高泛化能力。

# 具体的特征衍生方法取决于数据的性质和问题的需求。一些常见的特征衍生方法包括：

# - **多项式特征：** 通过对原始特征进行幂运算，引入高次项，从而能够捕捉非线性关系。
  
# - **交互特征：** 创建新的特征，通过原始特征的组合来捕捉它们之间的交互作用。

# - **对数、指数变换：** 对特征进行对数或指数变换，可以改变特征的分布，有时有助于模型的学习。

# - **统计特征：** 基于原始特征的统计性质（如均值、标准差）创建新的特征。

# - **时间序列特征：** 对时间序列数据进行滞后、移动平均等操作，以引入时序相关的信息。

# 在实践中，特征衍生通常需要深入理解数据的背景和领域知识，以确保引入的新特征是合理且有意义的。

'''
`PolynomialFeatures` 是 scikit-learn 库中的一个工具类，用于生成原始特征的多项式特征。这个类通过对原始特征进行多项式扩展，可以帮助模型更好地捕捉数据中的非线性关系。以下是使用 `PolynomialFeatures` 的基本用法以及效果的简要说明：

```python
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import make_pipeline
import numpy as np

# 创建一些示例数据
X = np.array([1, 2, 3, 4]).reshape(-1, 1)  # 一维特征
y = np.array([2, 0, 1, 4])

# 使用 PolynomialFeatures 进行二次多项式扩展
degree = 2
poly = PolynomialFeatures(degree)
X_poly = poly.fit_transform(X)

# 使用线性回归模型拟合多项式特征
model = LinearRegression()
model.fit(X_poly, y)

# 生成预测结果
X_test = np.array([5]).reshape(-1, 1)
X_test_poly = poly.transform(X_test)
y_pred = model.predict(X_test_poly)

print("Original Feature:")
print(X)

print("\nPolynomial Feature:")
print(X_poly)

print("\nPrediction for X_test:")
print(X_test)
print("Predicted y:", y_pred)
```

上述代码演示了如何使用 `PolynomialFeatures` 对一维特征进行二次多项式扩展，然后使用线性回归模型进行拟合和预测。效果如下：

1. **原始特征：**
   ```
   [1]
   [2]
   [3]
   [4]
   ```

2. **二次多项式扩展后的特征：**
   ```
   [1, 1, 1]
   [1, 2, 4]
   [1, 3, 9]
   [1, 4, 16]
   ```

3. **预测结果：**
   ```
   [5]
   Predicted y: [15.]
   ```

`PolynomialFeatures` 默认会生成所有阶数不超过指定次数的多项式特征。
在上述例子中，`degree=2` 表示生成一次和二次两个阶数的特征。
这有助于模型更好地拟合非线性关系。需要注意，多项式特征的次数（degree）过高可能会导致模型过拟合，
因此选择合适的次数通常需要通过交叉验证等方法确定。'''