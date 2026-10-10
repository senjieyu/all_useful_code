import pandas as pd
import numpy as np
import pymrmr
from typing import Union, List, Optional, Tuple

def mrmr_feature_selection(
    input_path: str,
    output_path: str,
    target_column: Optional[str] = None,
    feature_columns: Optional[Union[List[str], str, Tuple[slice, ...]]] = None,
    sheet_name: Union[str, int] = 0,
    method: str = 'MIQ',
    num_features: int = 40,
    drop_constants: bool = True,
    drop_non_numeric: bool = True,
    verbose: bool = True
) -> Union[List[str], None]:
    """
    基于你的工作代码优化的mRMR特征选择函数
    
    参数:
        input_path: 输入文件路径(支持Excel或CSV)
        output_path: 结果输出路径
        target_column: 目标列名(如果为None，假设第一列是目标列)
        feature_columns: 特征列选择方式，支持:
                        - None或'all': 使用所有其他列
                        - 列名列表: 如['feat1', 'feat2']
                        - 列名字符串: 单个列名如'feat1'
                        - iloc方式: 如(slice(2,10),)选择第2-9列
        sheet_name: Excel工作表名(仅Excel需要)
        method: mRMR方法('MIQ'或'MID')
        num_features: 要选择的特征数量
        drop_constants: 是否删除常数列
        drop_non_numeric: 是否删除非数值列
        verbose: 是否打印处理过程
    
    返回:
        选中的特征列表(失败返回None)
    """
    try:
        # 1. 数据读取
        if verbose:
            print(f"\n{'='*50}")
            print(f"开始读取数据: {input_path}")
            
        if input_path.endswith(('.xlsx', '.xls')):
            data = pd.read_excel(input_path, sheet_name=sheet_name)
        else:
            data = pd.read_csv(input_path, low_memory=False)
            
        if verbose:
            print(f"原始数据形状: {data.shape}")
            print("所有列:", list(data.columns))

        # 2. 处理特征列选择
        all_columns = data.columns.tolist()
        
        # 如果没有指定目标列，假设第一列是目标列
        if target_column is None:
            target_column = all_columns[0]
            if verbose:
                print(f"未指定目标列，使用第一列作为目标列: '{target_column}'")
        
        # 验证目标列存在
        if target_column not in data.columns:
            raise ValueError(f"目标列'{target_column}'不存在于数据中")
        
        if feature_columns is None or feature_columns == 'all':
            # 使用除目标列外的所有列
            feature_columns = [col for col in all_columns if col != target_column]
        elif isinstance(feature_columns, tuple) and len(feature_columns) == 1 and isinstance(feature_columns[0], slice):
            # iloc方式选择
            slice_obj = feature_columns[0]
            cols_idx = range(len(all_columns))[slice_obj]
            feature_columns = [all_columns[i] for i in cols_idx if all_columns[i] != target_column]
            
            if verbose:
                print(f"使用iloc方式选择的列索引: {list(cols_idx)}")
        elif isinstance(feature_columns, str):
            # 单个列名
            feature_columns = [feature_columns]
        elif isinstance(feature_columns, list):
            # 列名列表 - 保持不变
            pass
        else:
            raise ValueError(f"不支持的feature_columns类型: {type(feature_columns)}")
        
        # 验证特征列存在
        missing_cols = set(feature_columns) - set(all_columns)
        if missing_cols:
            raise ValueError(f"指定的特征列不存在: {missing_cols}")
            
        # 移除可能包含的目标列
        feature_columns = [col for col in feature_columns if col != target_column]
        
        if len(feature_columns) == 0:
            raise ValueError("没有可用的特征列")
        
        if verbose:
            print(f"\n目标列: '{target_column}'")
            print(f"初始特征列({len(feature_columns)}个): {feature_columns[:5]}...")  # 只显示前5个

        # 3. 根据pymrmr要求重新排列数据（目标列必须在第一列）
        selected_columns = [target_column] + feature_columns
        data = data[selected_columns]
        
        if verbose:
            print(f"重排后数据形状: {data.shape}")

        # 4. 删除常数列（参考你的代码）
        if drop_constants:
            constant_columns = data.columns[data.nunique() <= 1]  # 包括全NaN的情况
            if len(constant_columns) > 0:
                # 不删除目标列，即使它是常数
                constant_columns = [col for col in constant_columns if col != target_column]
                if constant_columns:
                    data = data.drop(columns=constant_columns)
                    if verbose:
                        print(f"删除常数列: {constant_columns}")

        # 5. 删除非数值列（参考你的代码）
        if drop_non_numeric:
            non_numeric_columns = data.select_dtypes(exclude=['number']).columns
            if len(non_numeric_columns) > 0:
                # 如果目标列是非数值的，需要特殊处理
                if target_column in non_numeric_columns:
                    if verbose:
                        print(f"目标列'{target_column}'是非数值类型，尝试转换...")
                    try:
                        # 尝试转换为数值
                        data[target_column] = pd.to_numeric(data[target_column], errors='coerce')
                        if data[target_column].isnull().all():
                            # 如果转换失败，使用标签编码
                            from sklearn.preprocessing import LabelEncoder
                            le = LabelEncoder()
                            data[target_column] = le.fit_transform(data[target_column].fillna('missing').astype(str))
                    except:
                        # 使用标签编码
                        from sklearn.preprocessing import LabelEncoder
                        le = LabelEncoder()
                        data[target_column] = le.fit_transform(data[target_column].fillna('missing').astype(str))
                
                # 删除其他非数值列
                non_numeric_to_drop = [col for col in non_numeric_columns if col != target_column]
                if non_numeric_to_drop:
                    data = data.drop(columns=non_numeric_to_drop)
                    if verbose:
                        print(f"删除非数值列: {non_numeric_to_drop}")

        if verbose:
            print(f"预处理后数据形状: {data.shape}")
            print("数据类型检查:")
            print(data.dtypes)

        # 6. 处理缺失值和无限值
        # 替换无限值
        data = data.replace([np.inf, -np.inf], np.nan)
        
        # 填充缺失值
        if data.isnull().sum().sum() > 0:
            if verbose:
                print("发现缺失值，正在处理...")
            # 数值列用中位数填充
            numeric_cols = data.select_dtypes(include=['number']).columns
            for col in numeric_cols:
                if data[col].isnull().sum() > 0:
                    fill_value = data[col].median()
                    if pd.isna(fill_value):
                        fill_value = 0
                    data[col] = data[col].fillna(fill_value)

        # 7. 确保数据类型正确（参考你的代码风格）
        if verbose:
            print("数据类型检查：")
            print("data 类型:", type(data))
            print("num_features 类型:", type(int(num_features)))

        # 确保num_features是整数且不超过特征数量
        num_features = int(num_features)
        max_features = data.shape[1] - 1  # 减去目标列
        if num_features > max_features:
            if verbose:
                print(f"警告: 请求的特征数({num_features})超过可用特征数({max_features})，调整为{max_features}")
            num_features = max_features

        # 8. 执行mRMR特征选择
        if verbose:
            print(f"\n执行mRMR特征选择(方法={method}, 特征数={num_features})...")
        
        try:
            selected_features = pymrmr.mRMR(data, method, num_features)
            if verbose:
                print("选中的特征：", selected_features)
        except Exception as e:
            print("执行 mRMR 时出错:", e)
            return None

        # 9. 保存结果（参考你的代码）
        if selected_features:
            # 创建包含目标列和选中特征的数据
            result_columns = [target_column] + selected_features
            # 从原始数据中提取（保持原始数据格式）
            if input_path.endswith(('.xlsx', '.xls')):
                original_data = pd.read_excel(input_path, sheet_name=sheet_name)
            else:
                original_data = pd.read_csv(input_path, low_memory=False)
            
            selected_data = original_data[result_columns]
            
            # 保存结果
            if output_path.endswith(('.xlsx', '.xls')):
                selected_data.to_excel(output_path, index=False)
            else:
                selected_data.to_csv(output_path, index=False)
                
            if verbose:
                print(f"选中的特征已保存到: {output_path}")
                print(f"最终数据形状: {selected_data.shape}")
                print('='*50)
        
        return selected_features
        
    except Exception as e:
        print(f"\n[错误] 处理过程中发生异常: {str(e)}")
        if verbose:
            import traceback
            traceback.print_exc()
        return None


# 使用示例（基于你的代码风格）
if __name__ == "__main__":
    # 示例1: 基本使用（类似你的原始代码）
    input_file = r"C:/Users/senji/Desktop/高光谱＋近红外/数据/近红外带标签合并结果.csv"
    output_file = r"C:/Users/senji/Desktop/高光谱＋近红外/数据/近红外带标签合并结果.csv"
    
    # print("示例1: 基本使用（假设第一列是目标列）")
    # selected_features = mrmr_feature_selection(
    #     input_path=input_file,
    #     output_path=output_file,
    #     target_column=None,  # 使用第一列作为目标列
    #     feature_columns=None,  # 使用所有其他列
    #     method='MIQ',
    #     num_features=40,
    #     verbose=True
    # )
    
    # # 示例2: 指定目标列和特征列范围
    # print("\n\n示例2: 指定目标列和特征列范围")
    # selected_features_2 = mrmr_feature_selection(
    #     input_path=input_file,
    #     output_path=output_file.replace(".xlsx", "_range.xlsx"),
    #     target_column="目标列名",  # 替换为实际的目标列名
    #     feature_columns=(slice(1, 100),),  # 选择第1-99列作为特征候选
    #     method='MIQ',
    #     num_features=50,
    #     verbose=True
    # )
    
    # 示例3: CSV文件处理
    print("\n\n示例3: CSV文件处理")
    csv_input = "C:/Users/senji/Desktop/高光谱＋近红外/数据/高光谱可用.csv"
    csv_output = "C:/Users/senji/Desktop/高光谱＋近红外/数据/高光谱100.csv"
    selected_features_csv = mrmr_feature_selection(
        input_path=csv_input,
        output_path=csv_output,
        target_column="part1",  # 指定目标列名
        feature_columns=(slice(4, 1000),),  # 选择第4-999列
        method='MIQ',
        num_features=100,
        verbose=True
    )