import subprocess
import os
import pandas as pd
from biom import load_table
import joblib

TAXONOMY = 'D:\\A_project\\code_system\\bioanalysis\\static\\Biom\\taxonomy_96.tsv'
BIOM = 'D:\\A_project\\code_system\\bioanalysis\\static\\Biom\\feature-table_96.biom'
RPART_MODEL = 'D:\\A_project\\code_system\\bioanalysis\\static\\model_parm\\Double_all_best_rpart_model.joblib'
RF_MODEL = 'D:\\A_project\\code_system\\bioanalysis\\static\\model_parm\\Double_all_best_rf_model.joblib'

def run_fastp(input_fastq, output_dir, log_file='media/analysis.log'):
    """
    运行 fastp 进行质控，并实时更新进度到指定文件
    :param input_fastq: 输入文件名（不包括 _1.fastq 和 _2.fastq 后缀）
    :param output_dir: 输出目录
    :param log_file: 日志文件路径
    :return: 如果成功运行返回 True, 否则返回 False
    """
    try:
        cmd = [
            "docker", "run", "--rm", "-v", f"{output_dir}:/data", "fastp-full",
            "-i", f"/data/{input_fastq}_1.fastq",
            "-I", f"/data/{input_fastq}_2.fastq",
            "-o", f"/data/{input_fastq}_clean_1.fastq",
            "-O", f"/data/{input_fastq}_clean_2.fastq",
            "--detect_adapter_for_pe",
            "--cut_front", "--cut_tail",
            "--length_required", "100",
            "-h", "/data/report.html",
            "-j", "/data/report.json",
            "-w", "4"
        ]
        with open(log_file, 'a', encoding='utf-8') as lf:
            process = subprocess.Popen(cmd, stdout=lf, stderr=lf)
            return_code = process.wait()  # 等待进程完成并获取返回码
            
            if return_code != 0:
                lf.write(f"Error: fastp process failed with return code {return_code}\n")
                return False
                
    except Exception as e:
        with open(log_file, 'a', encoding='utf-8') as lf:
            lf.write(f"Error running fastp: {e}\n")
        return False
    return True



def generate_manifest(data_dir, manifest_file):
    with open(manifest_file, 'w') as manifest:
        manifest.write("sample-id,absolute-filepath,direction\n")

        for file_name in os.listdir(data_dir):
            if file_name.endswith("_clean_1.fastq"):
                sample_name = file_name.replace("_clean_1.fastq", "")
                r2_file = f"{sample_name}_clean_2.fastq"

                manifest.write(f"{sample_name},/data/{file_name},forward\n")

                if os.path.exists(os.path.join(data_dir, r2_file)):
                    manifest.write(f"{sample_name},/data/{r2_file},reverse\n")
                else:
                    print(f"Warning: Missing {r2_file}, skipping reverse entry.")

    print(f"Manifest file generated: {manifest_file}")


def run_qiime2_analysis(data_dir, log_file='media/analysis.log'):
    try:
        model_dir = 'D:\\A_project\\code_system\\code_system\\models'
        model_file = model_dir + '\\' + 'silva-138-99-nb-classifier.qza'
        if not os.path.exists(model_file):
            with open(log_file, 'a', encoding='utf-8') as lf:
                lf.write(f"Error: Model file not found at {model_file}\n")
                lf.write("Please place silva-138-99-nb-classifier.qza in the models directory.\n")
            return False
        
        # 添加 QIIME2 分析开始日志
        with open(log_file, 'a', encoding='utf-8') as lf:
            lf.write("[INFO] Starting QIIME2 analysis...\n")
            lf.flush()
        
        cmd = [
            "docker", "run", "--rm", 
            "-v", f"{data_dir}:/data", 
            "-v", f"{model_dir}:/models",
            "quay.io/qiime2/amplicon:2025.4",
            "bash", "-c",
            "qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path /data/manifest.txt --output-path /data/paired-end-demux.qza --input-format PairedEndFastqManifestPhred33 && "
            "qiime demux summarize --i-data /data/paired-end-demux.qza --o-visualization /data/demux.qzv && "
            "qiime dada2 denoise-paired --i-demultiplexed-seqs /data/paired-end-demux.qza --p-trim-left-f 0 --p-trim-left-r 0 --p-trunc-len-f 0 --p-trunc-len-r 0 --o-table /data/table.qza "
            "--o-representative-sequences /data/rep-seqs.qza --o-denoising-stats /data/denoising-stats.qza && "
            "qiime metadata tabulate --m-input-file /data/denoising-stats.qza --o-visualization /data/denoising-stats.qzv && "
            "qiime feature-table summarize --i-table /data/table.qza --o-visualization /data/3_table_summary.qzv && "
            "qiime feature-table tabulate-seqs --i-data /data/rep-seqs.qza --o-visualization /data/rep-seqs.qzv && "
            "qiime phylogeny align-to-tree-mafft-fasttree --i-sequences /data/rep-seqs.qza --o-alignment /data/aligned-rep-seqs.qza "
            "--o-masked-alignment /data/masked-aligned-rep-seqs.qza --o-tree /data/unrooted-tree.qza --o-rooted-tree /data/rooted-tree.qza && "
            "qiime feature-classifier classify-sklearn --i-classifier /models/silva-138-99-nb-classifier.qza --i-reads /data/rep-seqs.qza --o-classification /data/taxonomy.qza && "
            "qiime metadata tabulate --m-input-file /data/taxonomy.qza --o-visualization /data/taxonomy.qzv && "
            "qiime tools export --input-path /data/table.qza --output-path /data/Biom && "
            "biom convert -i /data/Biom/feature-table.biom -o /data/Biom/OTU_table.tsv --to-tsv && "
            "qiime tools export --input-path /data/taxonomy.qza --output-path /data/Biom"
        ]
        
        with open(log_file, 'a', encoding='utf-8') as lf:
            lf.write("[INFO] Running QIIME2 commands...\n")
            lf.flush()
            
            # 使用实时输出而不是等待完成
            process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, bufsize=1, universal_newlines=True)
            
            # 实时读取输出并写入日志
            while True:
                output = process.stdout.readline()
                if output == '' and process.poll() is not None:
                    break
                if output:
                    lf.write(f"[QIIME2] {output.strip()}\n")
                    lf.flush()
            
            # 检查返回码
            return_code = process.poll()
            if return_code == 0:
                lf.write("[SUCCESS] QIIME2 analysis completed successfully.\n")
                lf.flush()
            else:
                lf.write(f"[ERROR] QIIME2 analysis failed with return code {return_code}\n")
                lf.flush()
                return False
                
    except subprocess.CalledProcessError as e:
        with open(log_file, 'a', encoding='utf-8') as lf:
            lf.write(f"[ERROR] Error running QIIME2 analysis: {e}\n")
            lf.flush()
        return False
    except Exception as e:
        with open(log_file, 'a', encoding='utf-8') as lf:
            lf.write(f"[ERROR] Unexpected error in QIIME2 analysis: {str(e)}\n")
            lf.flush()
        return False
    return True



def pred_model(data_dir):
    # 参考数据文件（固定路径，作为模板数据）
    taxonomy_file = TAXONOMY
    biom_file = BIOM

    if not os.path.exists(taxonomy_file) or not os.path.exists(biom_file):
        return False

    # 用户样本数据文件（来自传入的data_dir）
    taxonomy_file_one = data_dir + '\\Biom\\taxonomy.tsv'
    biom_file_one = data_dir + '\\Biom\\feature-table.biom'

    if not os.path.exists(taxonomy_file_one) or not os.path.exists(biom_file_one):
        return False
    
    # 读取taxonomy.tsv文件，第一行为表头
    taxonomy_df = pd.read_csv(taxonomy_file, sep='\t')
    taxonomy_df_one = pd.read_csv(taxonomy_file_one, sep='\t')

    # 加载biom文件
    table = load_table(biom_file)
    table_one = load_table(biom_file_one)

    # 提取计数矩阵，行是OTU/feature，列是样本
    counts = table.matrix_data.toarray()
    counts_one = table_one.matrix_data.toarray()

    # 行和列名称
    obs_ids = table.ids(axis='observation')  # OTU/feature ids
    sample_ids = table.ids(axis='sample')  # sample ids

    obs_ids_one = table_one.ids(axis='observation')
    sample_ids_one = table_one.ids(axis='sample')

    # 转成DataFrame: 行是OTU，列是样本
    counts_df = pd.DataFrame(counts, index=obs_ids, columns=sample_ids)
    counts_df_one = pd.DataFrame(counts_one, index=obs_ids_one, columns=sample_ids_one)

    amp_biom = counts_df.T
    amp_biom_one = counts_df_one.T

    # 创建一个列表用于存储每个参数的合并结果
    merged_columns = []

    # 按照每个参数将OTU列进行分组并相加-主体
    for parameter in taxonomy_df['Taxon'].unique():
        # 找出与当前参数对应的OTU（特征）
        relevant_features = taxonomy_df[taxonomy_df['Taxon'] == parameter]['Feature ID']

        # 获取这些特征在amp_biom中的列
        relevant_columns = amp_biom[relevant_features]

        # 将这些列的数值相加，并将结果存储在列表中
        merged_columns.append(relevant_columns.sum(axis=1))

    # 使用pd.concat一次性合并所有列，避免性能问题
    merged_df = pd.concat(merged_columns, axis=1)
    # 为每个参数列命名
    merged_df.columns = taxonomy_df['Taxon'].unique()
    dat = merged_df.div(merged_df.sum(axis=1), axis=0) * 100  # 按行归一化


    # 创建一个列表用于存储每个参数的合并结果--单个样本
    merged_columns_one = []

    # 按照每个参数将OTU列进行分组并相加-主体
    for parameter_one in taxonomy_df_one['Taxon'].unique():
        # 找出与当前参数对应的OTU（特征）
        relevant_features_one = taxonomy_df_one[taxonomy_df_one['Taxon'] == parameter_one]['Feature ID']

        # 获取这些特征在amp_biom中的列
        relevant_columns_one = amp_biom_one[relevant_features_one]

        # 将这些列的数值相加，并将结果存储在列表中
        merged_columns_one.append(relevant_columns_one.sum(axis=1))

    # 使用pd.concat一次性合并所有列，避免性能问题
    merged_df_one = pd.concat(merged_columns_one, axis=1)
    # 为每个参数列命名
    merged_df_one.columns = taxonomy_df_one['Taxon'].unique()
    dat_one = merged_df_one.div(merged_df_one.sum(axis=1), axis=0) * 100  # 按行归一化


    # 比较dat和dat_one, 将dat_one中缺少的列进行补零
    # 选择amp_biom的第一行
    dat_row = dat.iloc[0]

    # 创建一个空的DataFrame，用于保存结果
    dat_filled = pd.DataFrame(0.0, index=['sample'], columns=dat.columns)

    # 遍历dat_one的列，填充dat_filled中的对应列
    for column in dat_one.columns:
        if column in dat.columns:
            # 将dat_one中的数据填入dat_filled对应的列
            dat_filled[column].values[0] = dat_one[column].values[0]

    X = dat_filled

    # 读取模型
    rpart_model = joblib.load(RPART_MODEL)
    rf_model = joblib.load(RF_MODEL)

    # 测试模型表现（AUC评分）
    y_pred_rpart = rpart_model.predict_proba(X)
    y_pred_rf = rf_model.predict_proba(X)

    return y_pred_rpart, y_pred_rf