# %%
import re
import csv

def extract_gene_info_from_fasta(file_path):
    results = []  
    
    with open(file_path, 'r') as file:
        lines = file.readlines()  # 读取FASTA文件的所有行

    fasta_string = ''  # 存储当前基因的序列
    gene_id = None  # 存储当前基因的ID

    for line in lines:
        line = line.strip()  # 去除每行的前后空白字符
        if line.startswith('>'):
            # 如果遇到新的描述行（以'>'开头），先处理之前的序列
            if gene_id is not None:
                results.append(process_sequence(gene_id, fasta_string))  # 处理并保存之前的基因信息
            gene_id = line[1:].split(' ')[0]  # 提取基因ID（去除'>'和其他信息）
            fasta_string = ''  # 重置序列
        else:
            fasta_string += line  # 将当前行的序列拼接起来

    # 处理文件中的最后一个基因
    if gene_id is not None:
        results.append(process_sequence(gene_id, fasta_string))

    return results  # 返回提取的基因信息列表

def calculate_base_ratios(sequence):
    """
    计算序列中各个碱基(A, T, C, G)的比例。
    """
    length = len(sequence)  # 计算序列的总长度
    if length == 0:
        return {'A': 'NA', 'T': 'NA', 'C': 'NA', 'G': 'NA', 'A_ratio': 'NA', 'T_ratio': 'NA', 'C_ratio': 'NA', 'G_ratio': 'NA'}

    # 计算各个碱基的数量
    a_count = sequence.count('A') + sequence.count('a')
    t_count = sequence.count('T') + sequence.count('t')
    c_count = sequence.count('C') + sequence.count('c')
    g_count = sequence.count('G') + sequence.count('g')

    # 计算比例
    a_ratio = a_count / length
    t_ratio = t_count / length
    c_ratio = c_count / length
    g_ratio = g_count / length

    return {
        'A': a_count if a_count > 0 else 'NA',
        'T': t_count if t_count > 0 else 'NA',
        'C': c_count if c_count > 0 else 'NA',
        'G': g_count if g_count > 0 else 'NA',
        'A_ratio': a_ratio if a_ratio > 0 else 'NA',
        'T_ratio': t_ratio if t_ratio > 0 else 'NA',
        'C_ratio': c_ratio if c_ratio > 0 else 'NA',
        'G_ratio': g_ratio if g_ratio > 0 else 'NA'
    }

def process_sequence(gene_id, sequence):
    """
    计算5' UTR、CDS和3' UTR中碱基的比例。
    """
    # 查找5' UTR（前段小写字母）
    five_utr_match = re.match(r'^([a-z]+)', sequence)
    five_utr = five_utr_match.group(0) if five_utr_match else ''

    # 查找CDS（中段大写字母）
    cds_match = re.search(r'([A-Z]+)', sequence)
    cds = cds_match.group(0) if cds_match else ''

    # 查找3' UTR（后段小写字母）
    three_utr_match = re.search(r'([a-z]+)$', sequence)
    three_utr = three_utr_match.group(0) if three_utr_match else ''

    # 计算每个部分中的碱基比例
    five_utr_ratios = calculate_base_ratios(five_utr)
    cds_ratios = calculate_base_ratios(cds)
    three_utr_ratios = calculate_base_ratios(three_utr)

    return {
        'GeneID': gene_id,
        '5\'UTR_A': five_utr_ratios['A'],
        '5\'UTR_T': five_utr_ratios['T'],
        '5\'UTR_C': five_utr_ratios['C'],
        '5\'UTR_G': five_utr_ratios['G'],
        '5\'UTR_A_ratio': five_utr_ratios['A_ratio'],
        '5\'UTR_T_ratio': five_utr_ratios['T_ratio'],
        '5\'UTR_C_ratio': five_utr_ratios['C_ratio'],
        '5\'UTR_G_ratio': five_utr_ratios['G_ratio'],
        'CDS_A': cds_ratios['A'],
        'CDS_T': cds_ratios['T'],
        'CDS_C': cds_ratios['C'],
        'CDS_G': cds_ratios['G'],
        'CDS_A_ratio': cds_ratios['A_ratio'],
        'CDS_T_ratio': cds_ratios['T_ratio'],
        'CDS_C_ratio': cds_ratios['C_ratio'],
        'CDS_G_ratio': cds_ratios['G_ratio'],
        '3\'UTR_A': three_utr_ratios['A'],
        '3\'UTR_T': three_utr_ratios['T'],
        '3\'UTR_C': three_utr_ratios['C'],
        '3\'UTR_G': three_utr_ratios['G'],
        '3\'UTR_A_ratio': three_utr_ratios['A_ratio'],
        '3\'UTR_T_ratio': three_utr_ratios['T_ratio'],
        '3\'UTR_C_ratio': three_utr_ratios['C_ratio'],
        '3\'UTR_G_ratio': three_utr_ratios['G_ratio']
    }

def save_results_to_csv(results, output_file):
    """
    将结果保存到CSV文件中。
    """
    # 定义CSV的列名
    keys = results[0].keys() if results else [
        'GeneID', '5\'UTR_A', '5\'UTR_T', '5\'UTR_C', '5\'UTR_G', 
        '5\'UTR_A_ratio', '5\'UTR_T_ratio', '5\'UTR_C_ratio', '5\'UTR_G_ratio', 
        'CDS_A', 'CDS_T', 'CDS_C', 'CDS_G', 
        'CDS_A_ratio', 'CDS_T_ratio', 'CDS_C_ratio', 'CDS_G_ratio',
        '3\'UTR_A', '3\'UTR_T', '3\'UTR_C', '3\'UTR_G', 
        '3\'UTR_A_ratio', '3\'UTR_T_ratio', '3\'UTR_C_ratio', '3\'UTR_G_ratio'
    ]
    with open(output_file, 'w', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=keys)  # 创建CSV写入器
        writer.writeheader()  # 写入CSV头部
        writer.writerows(results)  # 写入所有行的数据

def main():
    fasta_file_path = '/Users/yding/Downloads/9311.IGDBv1.Allset.trans.fasta'  # FASTA文件路径
    csv_file_path = '/Users/yding/Downloads/genes_info.csv'  # 输出CSV文件路径

    results = extract_gene_info_from_fasta(fasta_file_path)  # 提取基因信息
    save_results_to_csv(results, csv_file_path)  # 保存结果到CSV文件
    print(f"Results saved to {csv_file_path}")  # 输出保存结果的提示信息

if __name__ == '__main__':
    main() 



