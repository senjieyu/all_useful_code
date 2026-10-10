
# 求导加平滑
#卢天代码
import pandas as pd
import numpy as np
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt

def savgol(
    spectra, # DataFrame, np.array
    window_length: int=7, 
    polyorder: int=2, 
    deriv: int=2
):
    spectra, _ = check_df(spectra)
    spectra = _savgol(spectra, window_length, polyorder, deriv)
    if _ is not None:
        spectra = _[0](spectra, **_[1])
    return spectra

def _savgol(
    spectra: np.ndarray, 
    window_length: int=7, 
    polyorder: int=2, 
    deriv: int=2
):
    return savgol_filter(spectra, window_length=window_length, polyorder=polyorder, deriv=deriv)

def check_df(data):
    if isinstance(data, pd.DataFrame):
        columns = data.columns
        index = data.index
        values = data.values
        return values, (pd.DataFrame, {'columns': columns, 'index': index})
    elif isinstance(data, np.ndarray):
        return data, None
    elif isinstance(data, pd.Series):
        index = data.index
        name = data.name
        values = data.values
        return values, (pd.Series, {'index': index, 'name': name})
    elif data is None:
        return None, None


def smooth_wavelengths(input_file_path,output_file):
    """
    读取CSV文件，对每一列的波长数据进行S-G平滑处理，并将处理后的数据保存为CSV文件
    :param input_file_path: 输入CSV文件路径
    :param output_file: 输出CSV文件路径
    """
    # 读取CSV文件
    df = pd.read_csv(input_file_path)
    
    # 获取波长列
    wavelength_columns = df.columns[2:]  # 假设波长数据从第三列开始
    
    # 取前两列数据
    df1 = df.iloc[:, :2]
    
    # 对波长数据进行S-G平滑处理
    df2 = savgol(df[wavelength_columns])
    
    # 合并处理后的数据与前两列数据
    smoothed_df = pd.concat([df1, df2], axis=1)
    
    # 将处理后的数据保存为CSV文件
    smoothed_df.to_csv(output_file)


input_file = r"C:/Users/senji/Desktop/梅菜扣肉/删除15异常值.csv"  # 输入CSV文件名
output_file = r"C:/Users/senji/Desktop/梅菜扣肉/一阶导平滑.csv"  # 输出CSV文件名
smooth_wavelengths(input_file, output_file)






# =============================================================================
#  S-G平滑的另一种方式（两种结果）
# import pandas as pd
# from scipy.signal import savgol_filter
# 
# def sg_smooth(data, window_size=5, poly_order=2):
#     """
#     使用S-G平滑处理给定的数据
#     :param data: 待处理的数据
#     :param window_size: 平滑窗口大小
#     :param poly_order: 多项式拟合阶数
#     :return: 平滑处理后的数据
#     """
#     return savgol_filter(data, window_size, poly_order)
# 
# def smooth_wavelengths(input_file, output_file):
#     """
#     读取CSV文件，对每一列的波长数据进行S-G平滑处理，并将处理后的数据保存为CSV文件
#     :param input_file: 输入CSV文件名
#     :param output_file: 输出CSV文件名
#     """
#     # 读取CSV文件
#     df = pd.read_csv(input_file)
#     
#     # 获取波长列
#     wavelength_columns = df.columns[2:]  # 假设波长数据从第三列开始
#     
#     # 取前两列数据
#     df1 = df.iloc[:, :2]
#     
#     # 对波长数据进行S-G平滑处理
#     df2 = df[wavelength_columns].apply(sg_smooth, window_size=5, poly_order=2)
#     
#     # 合并处理后的数据与前两列数据
#     smoothed_df = pd.concat([df1, df2], axis=1)
#     
#     # 将处理后的数据保存为CSV文件
#     smoothed_df.to_csv(output_file, index=False)
# 
# # 使用示例
# input_file = r"C:/Users/senji/Desktop/烟草仅部位/部位二阶导函数.csv"  # 输入CSV文件名
# output_file = r"C:/Users/senji/Desktop/烟草仅部位/二阶导平滑.csv"  # 输出CSV文件名
# smooth_wavelengths(input_file, output_file)
# 
# =============================================================================




















