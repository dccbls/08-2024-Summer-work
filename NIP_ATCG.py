# %%
import csv
from Bio import SeqIO

def calculate_base_ratios(sequence):
    """
    计算序列中各个碱基(A, T, C, G)的比例。
    """
    length = len(sequence) 
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

def process_sequence(gene_id, sequence, cds_start, cds_end):
    """
    处理每个基因的序列,计算5' UTR、CDS和3' UTR中碱基的比例。
    """
    # 注意SeqIO提供的序列是0-based的
    cds_start -= 1
    cds_end -= 1

    # 提取5'UTR、CDS 和 3'UTR
    five_utr = sequence[:cds_start]
    cds = sequence[cds_start:cds_end + 1]
    three_utr = sequence[cds_end + 1:]

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

def process_fasta_file(input_fasta, output_csv):
    """
    处理FASTA文件，计算每个基因的5'UTR、CDS和3'UTR区域中的碱基比例，并将结果保存到CSV文件中。
    """
    results = []

    for record in SeqIO.parse(input_fasta, "fasta"):
        header = record.description
        sequence = str(record.seq)

        if "CDS=" in header:
            # 提取CDS区域的起始和终止位置
            cds_range = header.split("CDS=")[-1]
            cds_start, cds_end = map(int, cds_range.split('-'))
            gene_id = header.split()[0]

            # 处理每个基因的序列
            result = process_sequence(gene_id, sequence, cds_start, cds_end)
            results.append(result)

    # 保存结果到CSV文件
    save_results_to_csv(results, output_csv)
    print(f"Results saved to {output_csv}")

input_fasta = "/Users/yding/Downloads/MSU.IGDBv1.Allset.mrna.fasta"  # 输入的FASTA文件路径
output_csv = "/Users/yding/NIP_ATCG.csv"         # 输出的CSV文件路径
process_fasta_file(input_fasta, output_csv)



