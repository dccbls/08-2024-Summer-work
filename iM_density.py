# %%
import csv
import pandas as pd
from collections import defaultdict

def process_data(file_path):
    id_count = defaultdict(int)  # 初始化计数器
    with open(file_path, 'r', encoding='utf-8') as file:
        # 读取文件的每一行
        for line in file:
            # 移除行尾的换行符并以制表符分割
            parts = line.strip().split('\t')
            # 使用前两个部分形成完整的ID，假设至少有两部分
            if len(parts) >= 2:
                original_id = parts[0] #+ parts[3]  # 合并第一部分和第二部分
                original_id = original_id.strip()  # 确保移除首尾的空白符
                id_count[original_id] += 1  # 对每个完整的ID进行计数

    # 将计数结果转为DataFrame
    id_list = list(id_count.keys())
    count_list = list(id_count.values())
    df = pd.DataFrame({
        'ID': id_list,
        '数量': count_list
    })

    return df

# 示例使用
file_path = '/Users/yding/Desktop/OS_conformation.txt'  # 替换为你的文件路径
df = process_data(file_path)

# 保存为Excel文件
output_file = '/Users/yding/Desktop/OOS_count.xlsx'  # 输出文件名
df.to_excel(output_file, index=False)

print(f"结果已保存为 {output_file}")


# %%
import pandas as pd

def parse_fasta(file_path):
    # 存储ID和序列长度的列表
    data = []
    
    # 读取FASTA文件
    with open(file_path, 'r') as f:
        seq_id = None
        seq = ""
        
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                # 如果已经有ID，保存上一个序列的信息
                if seq_id is not None:
                    data.append([seq_id, len(seq)])
                # 提取新的ID
                seq_id = line[1:]
                seq = ""
            else:
                # 拼接序列
                seq += line
        
        # 保存最后一个序列的信息
        if seq_id is not None:
            data.append([seq_id, len(seq)])
    
    # 转换为DataFrame
    df = pd.DataFrame(data, columns=["ID", "Length"])
    
    return df

# 示例使用
fasta_file = '/Users/yding/Documents/iwgsc_refseqv2.1_annotation_200916_HC_mrna.fasta'  # 替换为你的FASTA文件路径
output_excel = '/Users/yding/Desktop/mRNA_length.xlsx'

# 解析FASTA文件并生成Excel
df = parse_fasta(fasta_file)
df.to_excel(output_excel, index=False)

print(f"Excel表格已生成: {output_excel}")


# %%
import pandas as pd

# 读取Excel文件中的Sheet1
file_path = '/Users/yding/Desktop/CDS_length.xlsx'  # 替换为你的Excel文件路径
df = pd.read_excel(file_path, sheet_name='Sheet1')  # 只读取Sheet1

# 定义一个函数，根据ID提取组别（A、B、D等）
def get_group(id_value):
    if isinstance(id_value, str):
        # 提取ID中的组别字符，第8个字符表示组别
        return id_value[8]
    return None

# 应用函数，创建一个新的列来保存组别
df['Group'] = df['ID'].apply(get_group)

# 将数据按照组别分类，并保存到一个新的Excel文件中
output_file_path = '/Users/yding/Desktop/CDS_length_split.xlsx'
with pd.ExcelWriter(output_file_path, engine='openpyxl') as writer:
    for group_name, group_data in df.groupby('Group'):
        # 将每一组的数据保存到不同的工作表中
        group_data.to_excel(writer, sheet_name=f'Group_{group_name}', index=False)

print(f"分组结果已保存到: {output_file_path}")



