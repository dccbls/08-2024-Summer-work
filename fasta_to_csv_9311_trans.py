# %%
#import requirements
import re #导入正则表达式
import csv

def extract_gene_info_from_fasta(file_path):
    results = [] #用于储存所有的基因信息
    """
    从fasta文件中提取基因ID和各部分长度信息
    """
    with open(file_path, 'r') as file:
        lines = file.readlines() #读取fasta文件的所有行

    fasta_string = '' #储存当前基因序列
    gene_id = None #储存基因ID

    for line in lines:
        line = line.strip() #去除每行前后的空白字符
        if line.startswith('>'): #如果遇到‘>‘字符，先处理之前的序列
            if gene_id is not None:
                results.append(process_sequence(gene_id, fasta_string)) #处理并保存之前的基因信息
            gene_id = line[1:].split(' ')[0]  #提取基因ID
            fasta_string = '' #重置序列
        else:
            fasta_string += line #将当前行的序列拼接起来

    if gene_id is not None:
        results.append(process_sequence(gene_id, fasta_string))

    return results #返回提取的基因信息列表

def process_sequence(gene_id, sequence):
    """
    处理每个基因的序列,计算5'UTR,CDS,3'UTR的长度
    """
    five_utr_match = re.match(r'^([a-z]+)', sequence)
    five_utr_length = len(five_utr_match.group(0)) if five_utr_match else 0

    cds_match = re.search(r'([A-Z]+)', sequence)
    cds_length = len(cds_match.group(0)) if cds_match else 0

    three_utr_match = re.search(r'([a-z]+)$', sequence)
    three_utr_length = len(three_utr_match.group(0)) if three_utr_match else 0

    return {
        'GeneID': gene_id,
        '5\'UTR': five_utr_length,
        'CDS': cds_length,
        '3\'UTR': three_utr_length
    }

def save_results_to_csv(results, output_file):
    keys = results[0].keys() if results else ['GeneID', '5\'UTR', 'CDS', '3\'UTR']
    with open(output_file, 'w', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=keys)
        writer.writeheader()
        writer.writerows(results)

def main():
    fasta_file_path = '/Users/yding/Downloads/9311.IGDBv1.Allset.trans.fasta'
    csv_file_path = '/Users/yding/Downloads/genes_info.csv'

    results = extract_gene_info_from_fasta(fasta_file_path)
    save_results_to_csv(results, csv_file_path)
    print(f"Results saved to {csv_file_path}")

if __name__ == '__main__':
    main()



