import docker
import pandas as pd
import re
import logging
import argparse
import multiprocessing
from bs4 import BeautifulSoup
import shutil
import json
import os
import sys

from Fastq_QC import Fastq_QC
# %%
def run_docker(image, volumes, cmd):
    client = docker.from_env()
    return client.containers.run(image, cmd, volumes=volumes, remove=True, stdout=True, stderr=True)

def fq2asv(fq_list, outdir, samplename='test'):
    image_name='wnse/qiime2:20220718'
    cmd = f'python /home/qiime2/bin/Amplicon_fq2asv.py -i {" ".join(fq_list)} -o {outdir} -n {samplename}'
    vols = set([os.path.split(p)[0] for p in fq_list +[outdir]])
    vols = [f'{p}:{p}' for p in vols]    
    try:
        docker_out = run_docker(image_name, vols, cmd)
        logging.info(f'docker out: {docker_out}')
    except Exception as e:
        logging.error(e)
    outfile = os.path.join(outdir, 'Amplicon_fq2asv.json')
    if os.path.isfile(outfile):
        return outfile
    return None

def amplicon_merge(dir_list, outdir):
    image_name='wnse/qiime2:20220718'
    cmd = f'python /home/qiime2/bin/Amplicon_merge_with_diversity.py -i {" ".join(dir_list)} -o {outdir}'
    vols = set([os.path.split(p)[0] for p in dir_list + [outdir]])
    vols = [f'{p}:{p}' for p in vols]    
    try:
        docker_out = run_docker(image_name, vols, cmd)
        logging.info(f'docker out: {docker_out}')
    except Exception as e:
        logging.error(e)
    outfile = os.path.join(outdir, 'Amplicon_merge_with_diversity.json')
    if os.path.isfile(outfile):
        return outfile
    return None

def amplicon_diff(dir_list, metadata, outdir):
    image_name='wnse/qiime2:20220718'
    cmd = f'python /home/qiime2/bin/Amplicon_merge_with_diff.py -i {" ".join(dir_list)} -m {metadata} -o {outdir}'
    vols = set([os.path.split(p)[0] for p in dir_list + [metadata, outdir]])
    vols = [f'{p}:{p}' for p in vols]    
    try:
        docker_out = run_docker(image_name, vols, cmd)
        logging.info(f'docker out: {docker_out}')
    except Exception as e:
        logging.error(e)
    outfile = os.path.join(outdir, 'Amplicon_merge_with_diff.json')
    if os.path.isfile(outfile):
        return outfile
    return None

def parse_html_json(outdir):
    outdict = {}
    outdict['alpha_rarefaction_html_file'] = ''
    outdict['taxa_bar_html_file'] = ''

    alpha_rare_dir = os.path.join(outdir, 'alpha_rarefaction')
    if os.path.isdir(alpha_rare_dir):
        html_file = os.path.join(alpha_rare_dir, 'index.html')
        bak_file = os.path.join(alpha_rare_dir, 'index_bak.html')
        shutil.copy(html_file, bak_file)
        if os.path.isfile(html_file):
            try:
                with open(html_file,'rt') as h:
                    soup = BeautifulSoup(h, 'html.parser')
                soup.title.string = 'alpha_rarefaction'
                soup.find_all(id='q2templatesheader')[0].clear()
                with open(html_file,'w') as h:
                    h.write(soup.prettify(formatter='html'))
            except Exception as e:
                logging.error(f'parse alpha_rare_html {e}')
            outdict['alpha_rarefaction_html_file'] = alpha_rare_dir

    taxa_bar_dir = os.path.join(outdir, 'taxa_bar_plots')
    if os.path.isdir(taxa_bar_dir):
        html_file = os.path.join(taxa_bar_dir, 'index.html')
        bak_file = os.path.join(taxa_bar_dir, 'index_bak.html')
        shutil.copy(html_file, bak_file)
        if os.path.isfile(html_file):
            try:
                with open(html_file,'rt') as h:
                    soup = BeautifulSoup(h, 'html.parser')
                soup.title.string = 'barplot'
                soup.find_all(id='q2templatesheader')[0].clear()
                with open(html_file,'w') as h:
                    h.write(soup.prettify(formatter='html'))
            except Exception as e:
                logging.error(f'parse taxa_bar_html {e}')
            outdict['taxa_bar_html_file'] = taxa_bar_dir 

    return outdict

def parse_ancom_dict(meta_file, json_file):
    json_raw = {}
    outdict = []
    try:
        with open(json_file, 'rt') as H:
            json_raw = json.load(H)
        if 'ancom_res' in json_raw.keys():
            ancom_raw = json_raw['ancom_res']
    except Exception as e:
        logging.error(f'parse_ancom_dict {e}')

    df_abund = pd.DataFrame()
    try:
        taxa_infos = json_raw['taxa_info']
        for taxa_info in taxa_infos:
            # if taxa_info['level'] == 6:
            df_tmp = pd.DataFrame(taxa_info['info'])[taxa_info['legend']+['index']].set_index('index').T
            df_tmp['taxa_level'] = taxa_info['level']
            df_abund = pd.concat([df_abund, df_tmp])

        # if 'merged_tab_csv' in json_raw.keys():
        #     df = pd.read_csv(json_raw['merged_tab_csv'])
        #     df_abund = df.drop(['Sequence','Confidence'], axis=1).groupby('Taxon').sum()
        #     # df_abund = df_abund/df_abund.sum()
        #     df_abund.index = df_abund.index.str.replace(' ','')
        # else:
        #     logging.info(f' merged_tab_csv not in {json_file}')
    except Exception as e:
        logging.error(f'parse_tax_csv {e}')

    df_meta = pd.DataFrame()
    try:
        df_meta = pd.read_csv(meta_file, index_col=0, sep='\t')
    except Exception as e:
        logging.error(f'parse_meta_csv {e}')

    # try:
    for tag_group in ancom_raw:
        tmp_out = {}
        tag = ''
        raw_ancom_data = []
        if 'group' in tag_group.keys():
            tag = tag_group['group']
        if 'data' in tag_group.keys():
            raw_ancom_data = tag_group['data']
            df_data = pd.DataFrame(raw_ancom_data).set_index('id')

        tmp_out['group'] = tag
        tmp_out['data'] = raw_ancom_data
        tmp_out['ancom'] = []
        
        if ('ancom' in tag_group.keys()) and tag_group['ancom']:
            ancom_reject = pd.DataFrame(tag_group['ancom']).set_index('id')
            ancom_reject['clr'] = df_data.loc[ancom_reject.index,'clr']
            if tag in df_meta.columns:
                for tax in ancom_reject.index:
                    ancom_out = {}
                    ancom_out['id'] = tax
                    ancom_out['W'] = ancom_reject.loc[tax,'W']
                    ancom_out['clr'] = ancom_reject.loc[tax,'clr']
                    ancom_out['value'] = {}
                    ancom_out['value_dict'] = {}

                    level_len = len(tax.split(';'))
                    df_tmp_tax = df_abund[df_abund['taxa_level']==level_len].drop('taxa_level',axis=1).copy()
                    df_tmp_tax = df_tmp_tax/df_tmp_tax.sum()
                    # print(df_tmp_tax.to_dict())

                    for group in df_meta[tag].unique():
                        sample = df_meta[df_meta[tag]==group].index.to_list()
                        ancom_out['value'][group] = df_tmp_tax.loc[tax, sample].to_list()
                        ancom_out['value_dict'][group] = df_tmp_tax.loc[tax, sample].to_dict()
                tmp_out['ancom'].append(ancom_out)
            else:
                logging.info(f'{tag} not in {meta_file}')
        outdict.append(tmp_out)
    # except Exception as e:
        # logging.error(f'parse_ancom_json {e}')
    return outdict


def inchlib_json(data_csv, meta_csv, out_json):
    import inchlib_clust
    c = inchlib_clust.Cluster()
    c.read_csv(filename=data_csv, delimiter=",", header=True, missing_value=False, datatype="numeric")
    c.cluster_data(row_distance="euclidean", row_linkage="single", axis="both", column_distance="euclidean", column_linkage="single")
    d = inchlib_clust.Dendrogram(c)
    # create the cluster heatmap representation and define whether you want to compress the data by defining the maximum number of heatmap rows, the resulted value of compressed (merged) rows and whether you want to write the features
    # d.create_cluster_heatmap(compress=int, compressed_value="median", write_data=bool)
    d.create_cluster_heatmap()
    d.add_metadata_from_file(metadata_file=meta_csv, delimiter=",", header=True, metadata_compressed_value="frequency")
    # read column metadata file with specified delimiter, also specify whether there is a 'header' column
    # d.add_column_metadata_from_file(column_metadata_file="/path/to/file.csv", delimiter=",", header=bool)
    d.export_cluster_heatmap_as_json(out_json)


def get_inchlib_json(meta_file, json_file, outdir):
    json_raw = {}
    outdict = []
    beta_diversity = []
    # print(json_file)
    try:
        with open(json_file, 'rt') as H:
            json_raw = json.load(H)
        if 'beta_diversity' in json_raw.keys():
            beta_diversity = json_raw['beta_diversity']
    except Exception as e:
        logging.error(f'parse_beta_diversity_dict {e}')

    df_meta = pd.DataFrame()
    try:
        df_meta = pd.read_csv(meta_file, index_col=0, sep='\t')
    except Exception as e:
        logging.error(f'parse_meta_csv {e}')

    if beta_diversity:
        for group in df_meta.columns:
            tmp_out = {}
            tmp_out[group] = []
            for metric_lst in beta_diversity:
                metric = metric_lst['metric']
                tmp_out_metric = {}
                tmp_out_metric['metric'] = metric
                tmp_out_metric['matrix_clust_json'] = ''
                matrix_csv = metric_lst['distance_matrix_csv']
                
                meta_tmp_file = os.path.join(outdir,'tmp_matrix_meta.csv')
                df_meta[group].to_csv(meta_tmp_file)

                matrix_tmp_file = os.path.join(outdir,'tmp_matrix.csv')
                df_matrix = pd.read_csv(matrix_csv, index_col=0, sep='\t')
                df_matrix.index.name = df_meta.index.name
                df_matrix.to_csv(matrix_tmp_file)

                matrix_tmp_json = os.path.join(outdir,f'matrix_{metric}_{group}_clust.json')
                inchlib_json(matrix_tmp_file, meta_tmp_file, matrix_tmp_json)
                if os.path.isfile(matrix_tmp_json):
                    tmp_out_metric['matrix_clust_json'] = matrix_tmp_json
                tmp_out[group].append(tmp_out_metric)
            outdict.append(tmp_out)
    return outdict


if __name__ == '__main__':
    cpu_num = multiprocessing.cpu_count()
    if cpu_num > 4:
        cpu_num = cpu_num - 2
    bin_dir = os.path.split(os.path.realpath(__file__))[0]

    from mkdir import mkdir
    from write_json import write_json
    from post_status import post_url
    from post_status import post_pid
    from post_status import write_status
    
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', nargs='+', required=True, help='fastq file for assemble')
    parser.add_argument('-o', '--outdir', required=True, help='out dir for outputfiles')
    parser.add_argument('-m', '--metadata', default=None, help='meta csv for analysis (absolute path)')
    parser.add_argument('-n', '--name', default='test', help='sample name')
    parser.add_argument('-t', '--threads', default=cpu_num, help='threads for fastqc')
    parser.add_argument('-debug', '--debug', action='store_true')
    parser.add_argument('-type', '--type', default='fq2asv', choices=['fq2asv', 'merge'], help='choose function')
    parser.add_argument('-tag', '--tag', default='AsvAmplicon', help='output json file name')
    parser.add_argument('-tID', '--taskID', default='', help='task ID for report status')
    args = parser.parse_args()
    outdir = args.outdir
    logfile = os.path.join(outdir, 'log')
    mkdir(outdir)
    logging.basicConfig(level=logging.INFO, filename=logfile, format='%(asctime)s %(levelname)s %(message)s',datefmt='%Y-%m-%d %H:%M:%S')
    input_list = args.input
    sample_name = args.name
    metadata = args.metadata
    tag = args.tag

    info_dict = {}
    taskID = args.taskID
    post_pid(taskID)
    status_report = os.path.join(outdir, 'status_report.txt')

    if args.type == 'merge':
        step = tag
        logging.info(step)
        if not metadata:
            sys.exit('metadata should be given if run merge')
        if not os.path.isfile(metadata):
            sys.exit(f'metadata file not exists: {metadata}')
        try:
            s = f'{step}\tR\t'
            write_status(status_report, s)
            try:
                merge_json_file = amplicon_diff(input_list, metadata, outdir)
            except Exception as e:
                logging.error(e)

            # try:
            # shutil.copy(merge_json_file, os.path.join(outdir, f'{step}.json'))
            # os.symlink(merge_json_file, os.path.join(outdir, f'{step}.json'))
            # with open(os.path.join(outdir, f'{step}.json'),'rt') as H:
                # info_dict.update({'ASV_info':json.load(H)})
            with open(merge_json_file, 'rt') as H:
                tmp_out_dict = json.load(H)
            try:
                tmp_out_dict['html_files'] = parse_html_json(outdir)
                tmp_out_dict['ancom_res'] = parse_ancom_dict(metadata, merge_json_file)
                tmp_out_dict['beta_matrix_json'] = get_inchlib_json(metadata, merge_json_file, outdir)
            except Exception as e:
                logging.error(f'process_html_json {e}')
            with open(os.path.join(outdir, f'{step}.json'), 'w') as H:
                json.dump(tmp_out_dict, H, indent=2, ensure_ascii=False)

            # except Exception as e:
                # logging.error(f'symlink {merge_json_file} {e}')

            s = f'{step}\tD\t'
        except Exception as e:
            logging.error(e)
            s = f'{step}\tE\t'
        try:
            write_status(status_report, s)
            post_url(taskID, step)
        except Exception as e:
            logging.error(f'{step} status {e}')
        try:
            post_url(taskID, '2', 'http://localhost/task/getMergeStatus/')
        except Exception as e:
            logging.error(f'post_url getTaskRunningStatus {e}')

    if args.type == 'fq2asv':
        fq_list = input_list
        step = 'FastqcAmplicon'
        logging.info(step)
        try:
            s = f'{step}\tR\t'
            write_status(status_report, s)
            fq_cor_list, outdict = Fastq_QC(fq_list, outdir, threads=args.threads)
            if fq_cor_list:
                fq_list = fq_cor_list
                info_dict.update(outdict)
            with open(os.path.join(outdir, f'{step}.json'), 'w') as H:
                json.dump(outdict, H, indent=2)
            s = f'{step}\tD\t'
        except Exception as e:
            logging.error(e)
            s = f'{step}\tE\t'
        try:
            write_status(status_report, s)
            post_url(taskID, step)
        except Exception as e:
            logging.error(f'{step} status {e}')


        step = tag
        logging.info(step)
        try:
            s = f'{step}\tR\t'
            write_status(status_report, s)
            fq2asv_json_file = fq2asv(fq_list, outdir, samplename=args.name)
            if fq_cor_list:
                fq_list = fq_cor_list
                info_dict.update(outdict)

            try:
                shutil.copy(fq2asv_json_file, os.path.join(outdir, f'{step}.json'))
                # with open(os.path.join(outdir, f'{step}.json'),'rt') as H:
                    # info_dict.update({'ASV_info':json.load(H)})
            except Exception as e:
                logging.error(f'rename {fq2asv_json_file} {e}')

            s = f'{step}\tD\t'
        except Exception as e:
            logging.error(e)
            s = f'{step}\tE\t'
        try:
            write_status(status_report, s)
            post_url(taskID, step)
        except Exception as e:
            logging.error(f'{step} status {e}')

        try:
            post_url(taskID, '2', 'http://localhost/task/getTaskRunningStatus/')
        except Exception as e:
            logging.error(f'post_url getTaskRunningStatus {e}')

    # json_out = write_json(info_dict, outdir=outdir)
    # if not json_out:
        # logging.info(f'write json failed {info_dict}')
    if not args.debug:
        tmp_dir = os.path.join(outdir, 'tmp')
        try:
            if os.path.isdir(tmp_dir):
                shutil.rmtree(tmp_dir)
        except Exception as e:
            logging.error(e)

