# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.13.7
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %%
import pandas as pd
import json
import os
import logging
import argparse
import json
import re

# %%


def export_alpha_sig(vector_artifact, metadata_md):
    import qiime2.plugins.diversity.actions as diversity_actions
    viz, = diversity_actions.alpha_group_significance(
        alpha_diversity=vector_artifact,
        metadata=metadata_md,
    )
    return viz
    

def get_alpha_vectors_sig(alpha_index_csv, metadata_csv, outdir):
    from qiime2 import Artifact
    from qiime2 import Metadata
    metadata_md = Metadata.load(metadata_csv)
    df_vectors = pd.read_csv(alpha_index_csv, index_col=0)
    metrics = ['faith_pd', 'observed_features', 'shannon_entropy' ,'pielou_evenness']
    for m in metrics:
        if m in df_vectors.columns:
            tmp_vector = Artifact.import_data('SampleData[AlphaDiversity]', df_vectors[m])
            tmp_viz = export_alpha_sig(tmp_vector, metadata_md)
            tmp_viz.export_data(os.path.join(outdir, m))
        else:
            logging.info(f'{m} not in {alpha_index_csv} columns')
        

def export_ancom_sig(table_artifact, metadata_categroy):
    import qiime2.plugins.composition.actions as composition_actions
    comp_table, = composition_actions.add_pseudocount(
        table=table_artifact,
    )
    ancom_viz, = composition_actions.ancom(
        table=comp_table,
        metadata=metadata_categroy,
    )
    return ancom_viz

def get_ancom_sig(table, metadata, outdir):  
    for group in metadata.columns:
        try:
            ancom_viz = export_ancom_sig(table, metadata.get_column(group))
            ancom_viz.export_data(os.path.join(outdir, f'ancom_{group}'))
        except Exception as e:
            logging.error(e)
            
def get_ancom_sig_taxa(table_qza, metadata_csv, tax_qza, outdir):
    from qiime2 import Artifact
    from qiime2 import Metadata
    import qiime2.plugins.taxa.actions as taxa_actions
    metadata_md = Metadata.load(metadata_csv)
    table = Artifact.load(table_qza)
    taxonomy = Artifact.load(tax_qza)
    # for i in range(7):
    tax_table, = taxa_actions.collapse(
        table=table,
        taxonomy=taxonomy,
        level=6,
    )
    get_ancom_sig(tax_table, metadata_md, outdir)
    


# %%
def get_difference(merge_info, metadata_csv, outdir):
    outdict = {}
    if not os.path.isfile(metadata_csv):
        logging.error(f'file not exists {metadata_csv}')
        return outdict

    try:
        table_qza = merge_info['merged_tab_qza']
        tax_qza = merge_info['merged_tax_qza']
        if os.path.isfile(table_qza) and os.path.isfile(tax_qza):
            get_ancom_sig_taxa(table_qza, metadata_csv, tax_qza, outdir)
        else:
            logging.info(f'file not exists {table_qza} or {tax_qza}')
    except Exception as e:
        logging.error(e)

    try:
        alpha_index_csv = merge_info['alpha_index_csv']
        if os.path.isfile(alpha_index_csv):
            get_alpha_vectors_sig(alpha_index_csv, metadata_csv, outdir)
        else:
            logging.info(f'file not exists {alpha_index_csv}')
    except Exception as e:
        logging.error(e)



# %%

def check_file_exists(file, file_type='file'):
    if file_type == 'file':
        check = os.path.isfile
    elif file_type == 'dir':
        check = os.path.isdir
    else:
        return None
    return check(file)


def get_difference_ancom_res(meta_list, outdir):
    out = []
    for meta in meta_list:
        outdict = {}
        outdict['group'] = meta
        data_csv = os.path.join(outdir, f'ancom_{meta}', 'data.tsv')
        if check_file_exists(data_csv):
            outdict['data_csv'] = data_csv
            df = pd.read_csv(data_csv, sep='\t').astype(str)
            outdict['data'] = [v for i,v in df.T.to_dict().items()]

        ancom_csv = os.path.join(outdir, f'ancom_{meta}', 'ancom.tsv')
        reject_true_list = []
        if check_file_exists(ancom_csv):
            outdict['ancom_csv'] = ancom_csv
            df = pd.read_csv(ancom_csv, sep='\t').astype(str)
            df.columns = ['id', 'W', 'Reject']
            df = df[df['Reject']=='True']
            outdict['ancom'] = [v for i,v in df.T.to_dict().items()]
            if not df.empty:
                reject_true_list = df['id'].to_list()

        abund_csv = os.path.join(outdir, f'ancom_{meta}', 'percent-abundances.tsv')
        if check_file_exists(abund_csv):
            outdict['abund_csv'] = abund_csv
            df = pd.read_csv(abund_csv, sep='\t', index_col=0, header=None).astype(str)
            outdict['abund_percent'] = df.loc['Percentile'].to_list()
            outdict['abund_group'] = df.loc['Group'].to_list()
            df = df[df.index.isin(reject_true_list)]
            outdict['abund'] = df.apply(lambda x: {'id':x.name, 'value':x.to_list()}, axis=1).to_list()
            # outdict['abund'] = [{'id':i, 'value':v } for i, v in tmp_dict.items()]
        out.append(outdict)
    return out

def get_difference_alpha_res(meta_list, outdir):
    out = []
    metrics = ['faith_pd','observed_features','pielou_evenness','shannon_entropy']
    for meta in meta_list:
        outdict = {}
        outdict['group'] = meta
        outdict['metrics'] = []
        for metric in metrics:
            metric_tmp_dict = {}
            metric_tmp_dict['metric'] = metric
            jsonp_file = os.path.join(outdir, metric, f'column-{meta}.jsonp')
            if check_file_exists(jsonp_file):
                jsonp_content = open(jsonp_file, 'rt').readline()
                metric_tmp_dict['info'] = {}
                for tmp_dict in json.loads('['+ re.match(r'.*?({.*})',jsonp_content).group(1) +']'):
                    for i,v in tmp_dict.items():
                        metric_tmp_dict['info'][i] = v

            pairwise_file = os.path.join(outdir, metric, f'kruskal-wallis-pairwise-{meta}.csv')
            if check_file_exists(pairwise_file):
                pairwise_content = pd.read_csv(pairwise_file)
                metric_tmp_dict['pairwise'] = [v for i, v in pairwise_content.to_dict(orient='index').items()]

            if metric_tmp_dict:
                outdict['metrics'].append(metric_tmp_dict)
        out.append(outdict)
    return out


def get_difference_res(metadata_csv, outdir):
    outdict = {}
    meta_list = pd.read_csv(metadata_csv, sep='\t', index_col=0).columns
    outdict.update({'ancom_res':get_difference_ancom_res(meta_list, outdir)})
    outdict.update({'alpha_res':get_difference_alpha_res(meta_list, outdir)})
    return outdict



# %%

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s',datefmt='%Y-%m-%d %H:%M:%S')
    parse = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parse.add_argument('-i', '--input', required=True, help='merge json for analysis (absolute path)')
    parse.add_argument('-m', '--metadata', required=True, help='meta csv for analysis (absolute path)')
    parse.add_argument('-o', '--outdir', required=True, help='out dir for output files')
    args = parse.parse_args()
    
    from write_json import write_json
    from mkdir import mkdir

    json_file = args.input
    metadata_csv = args.metadata
    outdir = args.outdir
    if os.path.isfile(json_file):
        with open(json_file, 'rt') as H:
            merge_info = json.load(H)
        if os.path.isfile(metadata_csv):
            mkdir(outdir)
            info_dict = {}
            try:
                get_difference(merge_info, metadata_csv, outdir)
                info_dict = get_difference_res(metadata_csv, outdir)
            except Exception as e:
                logging.error(e)
            json_out = write_json(info_dict, outdir=outdir)
            if not json_out:
                logging.info(f'write json failed')
                logging.info(f'{info_dict}')
        else:
            logging.error(f'file not exists {metadata_csv}')
    else:
        logging.error(f'file not exists {json_file}')
    