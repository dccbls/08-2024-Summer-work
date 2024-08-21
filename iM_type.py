# %%
import pandas as pd

# 读取Excel文件
file_path = '/Users/yding/Documents/CS_mRNA_conformation.xlsx'  # 替换为你的文件路径
df = pd.read_excel(file_path)

# 定义分类函数
def categorize(row):
    category = None
    if row[3] == 3:
        max_L = max(row[4], row[5], row[6])
        if 1 <= max_L <= 4:
            category = 'C3L1-4'
        elif 5 <= max_L <= 8:
            category = 'C3L5-8'
        elif 9 <= max_L <= 12:
            category = 'C3L9-12'
    elif row[3] == 4:
        category = 'C4L1-12'
    elif row[3] >= 5:
        category = 'C5L1-12'
    
    return category

# 应用分类函数
df['Category'] = df.apply(categorize, axis=1)

# 统计每个类别的计数
category_counts = df['Category'].value_counts()

# 将计数结果转换为DataFrame
category_counts_df = category_counts.reset_index()
category_counts_df.columns = ['Category', 'Count']

# 输出为新的Excel文件
output_file_path = '/Users/yding/Documents/mRNA_categorized_counts.xlsx'  # 输出文件路径
category_counts_df.to_excel(output_file_path, index=False)

print(f'分类计数结果已保存到 {output_file_path}')



