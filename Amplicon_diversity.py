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
import sys
import os
import json
import pandas as pd
import logging
import argparse
import re


# %%
def get_diversity_action(sample_meta, table, tree):
    import qiime2.plugins.diversity.actions as diversity_actions
    import biom
    tmp_table = table.view(biom.Table)
    max_depth = int(min(tmp_table.sum(axis='sample')))
    action_results = diversity_actions.core_metrics_phylogenetic(
        phylogeny=tree,
        table=table,
        sampling_depth=max_depth,
        metadata=sample_meta,
    )
    alpha_rarefaction_viz, = diversity_actions.alpha_rarefaction(
        phylogeny=tree,
        table=table,
        max_depth=max_depth,
        metadata=sample_meta,
    )
    return action_results, alpha_rarefaction_viz

# %%

def check_file_exists(file, file_type='file'):
    if file_type == 'file':
        check = os.path.isfile
    elif file_type == 'dir':
        check = os.path.isdir
    else:
        return None
    return check(file)

def get_alpha_rarefaction_results(meta_list, alpha_rarefaction_viz, outdir):
    try:
        alpha_rarefaction_viz.export_data(os.path.join(outdir, f'alpha_rarefaction'))
    except Exception as e:
        logging.error(f'get_alpha_rarefaction_results:{e}')
    out = []
    metrics = {'observed_features', 'shannon', 'faith_pd'}
    for meta in meta_list:
        outdict = {}
        outdict['group'] = meta
        outdict['metrics'] = []
        for metric in metrics:
            metric_tmp_dict = {}
            metric_tmp_dict['metric'] = metric
            metric_tmp_dict['info'] = {}
            jsonp_file = os.path.join(outdir, 'alpha_rarefaction', f'{metric}-{meta}.jsonp')
            if check_file_exists(jsonp_file):
                jsonp_content = open(jsonp_file, 'rt').readline()
                for i,v in json.loads(re.match(r'.*?({.*})',jsonp_content).group(1)).items():
                    metric_tmp_dict['info'][i] = v
            if metric_tmp_dict:
                outdict['metrics'].append(metric_tmp_dict)
        out.append(outdict)
    return out

# %%
def get_diversity_pcoa_results(action_result, outdir):
    from skbio.stats.distance import DistanceMatrix
    from skbio import OrdinationResults
    
    outdict = {}
    df_alpha_idx = pd.concat([
        action_result.rarefied_table.view(pd.DataFrame).sum(axis=1).rename('sampleing_depth'),
        action_result.faith_pd_vector.view(pd.Series),
        action_result.observed_features_vector.view(pd.Series),
        action_result.shannon_vector.view(pd.Series),
        action_result.evenness_vector.view(pd.Series)],
        axis=1
    )
    alpha_idx_file = os.path.join(outdir, f'alpha_index.csv')
    df_alpha_idx.to_csv(alpha_idx_file)
    outdict.update({f'alpha_index_csv': alpha_idx_file})
    
    metrics = ['jaccard', 'bray_curtis',
               'unweighted_unifrac', 'weighted_unifrac']
    outdict['beta_diversity'] = []
    for metric in metrics:
        try:
            tmp_dict = write_Ord_Res(action_result, metric, outdir)
            tmp_dict['metric'] = metric
            outdict['beta_diversity'].append(tmp_dict)
        except Exception as e:
            logging.error(e)
    
    return outdict


# %%
def write_Ord_Res(action_result, metric, outdir):
    from skbio import OrdinationResults
    from skbio.stats.distance import DistanceMatrix
    outdict = {}
    
    distance_matrix = f'{metric}_distance_matrix'
    pcoa_result = f'{metric}_pcoa_results'
    
    sample_file = os.path.join(outdir, f'diversity_{metric}_sample.csv')
    exp_file = os.path.join(outdir, f'diversity_{metric}_proportion_explained.csv')
    dis_file = os.path.join(outdir, f'diversity_{metric}_distance_matrix.csv')
    
    df_samples = getattr(action_result, pcoa_result).view(OrdinationResults).samples
    df_exp = getattr(action_result, pcoa_result).view(OrdinationResults).proportion_explained
    df_dis = getattr(action_result, distance_matrix).view(DistanceMatrix).to_data_frame()
    
    df_samples.to_csv(sample_file, sep='\t')
    outdict.update({f'sample_pcoa_csv': sample_file})
    outdict.update({'sample_pcoa':[{"sample": i, "info": v} for i, v in df_samples.to_dict(orient="index").items()]})

    df_exp.to_csv(exp_file, header=None, sep='\t')
    outdict.update({f'pcoa_exp_csv': exp_file})
    outdict.update({'pcoa_exp':df_exp.to_dict()})
    
    df_dis.to_csv(dis_file, sep='\t')
    outdict.update({f'distance_matrix_csv':dis_file})
    outdict.update({'distance_matrix':df_dis.to_dict()})
    
    return outdict


# %%
def get_alpha_rarefaction(table, tree, outdir):
    outdict = {}
    from q2_diversity._alpha._visualizer import _compute_rarefaction_data
    from skbio import TreeNode
    import biom
    # from q2_diversity._alpha._visualizer import _compute_summary
    # from q2_diversity._alpha._visualizer import _alpha_rarefaction_jsonp
    outdict['alpha_rarefaction'] = {}
    metrics = {'observed_features', 'shannon', 'faith_pd'}
    tmp_table = table.view(biom.Table)
    max_depth = max(tmp_table.sum(axis='sample'))
    min_depth = 1
    steps = 10
    iterations = 10
    tmp_tree = tree.view(TreeNode)
    div_data = _compute_rarefaction_data(tmp_table, min_depth, max_depth, steps, iterations, tmp_tree, metrics)
    for m, data in div_data.items():
        filename = os.path.join(outdir, f'alpha_rarefaction_{m}.csv')
        # jsonp_filename = f'{m}.jsonp'
        # n_df = _compute_summary(data, 'sample-id')
        # _alpha_rarefaction_jsonp(outdir, jsonp_filename, m, n_df, '')
        data.columns = [f'depth_{t[0]}_iter{t[1]}' for t in data.columns.values]
        data.to_csv(filename)
        outdict['alpha_rarefaction'][m] = filename
    return outdict


# %%
def get_diversity(table_qza, tree_qza, outdir, meta=None):
    from qiime2 import Artifact
    from qiime2 import Metadata
    
    table = Artifact.load(table_qza)
    tree = Artifact.load(tree_qza)
    if not meta:
        tmpmanifest = os.path.join(outdir, 'tmpmanifest')
        df_tmp = pd.Series(table.view(pd.DataFrame).index.to_list()).rename('sample-id')
        df_tmp = pd.concat([df_tmp, df_tmp.rename('sample')], axis=1)
        df_tmp.to_csv(tmpmanifest, sep='\t', index=False)
    else:
        tmpmanifest = meta
    meta_list = pd.read_csv(tmpmanifest, sep='\t', index_col=0).columns
    sample_meta = Metadata.load(tmpmanifest)
    
    outdict = {}
    try:
        action_result, alpha_rarefaction_viz = get_diversity_action(sample_meta, table, tree)
        outdict = get_diversity_pcoa_results(action_result, outdir)
        outdict.update({'alpha_rarefaction':get_alpha_rarefaction_results(meta_list, alpha_rarefaction_viz, outdir)})
    except Exception as e:
        logging.error(e)
    # try:
        # outdict.update(get_alpha_rarefaction(table, tree, outdir))
    # except Exception as e:
    #     logging.error(e)
    return outdict


# %%
if __name__ == '__main__':
    bin_dir = os.path.split(os.path.realpath(__file__))[0]
    # pub_path = os.path.join(bin_dir, '../pub/')
    # if os.path.isdir(pub_path):
    #     sys.path.append(pub_path)
    # else:
    #     raise(f'{pub_path} not exists')
    
    from write_json import write_json
    from mkdir import mkdir
    
    parse = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parse.add_argument('-t', '--table', required=True, help='merged table qza file')
    parse.add_argument('-e', '--tree', required=True, help='merged rooted tree qza file')
    parse.add_argument('-m', '--meta', default=None, help='metadata csv ')
    parse.add_argument('-o', '--outdir', required=True, help='out dir for output files')
    args = parse.parse_args()
    
    table_qza = args.table
    tree_qza = args.tree
    
    outdir = args.outdir
    mkdir(outdir)
    logfile = os.path.join(outdir, 'log')
    logging.basicConfig(level=logging.INFO, filename=logfile, format='%(asctime)s %(levelname)s %(message)s',datefmt='%Y-%m-%d %H:%M:%S')
    
    if args.meta:
        if os.path.isfile(args.meta):
            metadata = args.meta
        else:
            metadata = None
    else:
        metadata = None

    info_dict = get_diversity(table_qza, tree_qza, outdir, meta=metadata)
    json_out = write_json(info_dict, outdir=outdir)
    if not json_out:
        logging.info(f'write json failed')
        logging.info(f'{info_dict}')

# %%
