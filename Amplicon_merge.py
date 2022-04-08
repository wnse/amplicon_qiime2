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


# %%
def data2df(data):
    from qiime2 import Metadata
    return data.view(Metadata).to_dataframe()


# %%
def get_qza_list(json_list, qza_name='asv_seq_qza'):
    from qiime2 import Artifact
    qza_list = []
    for j in json_list:
        with open(j,'rt') as h:
            fq2asv = json.load(h)
            if qza_name  in fq2asv.keys():
                file = fq2asv[qza_name]
                if os.path.isfile(file):
                    try:
                        qza_list.append(Artifact.load(file))
                    except Exceptions as e:
                        logging.error(e)
                else:
                    logging.info(f'{file} not exists')
            else:
                logging.info(f'{qza_name} not in {j}')
    return qza_list

def merge_qza(qza_list, outfile, para='seq'):
    from qiime2.plugins.feature_table.methods import merge
    from qiime2.plugins.feature_table.methods import merge_seqs
    from qiime2.plugins.feature_table.methods import merge_taxa
    from qiime2 import Metadata
    
    merge_para = None
    if para == 'tab':
        merge_para = merge
        merged = merge_para(qza_list)
        merged.merged_table.save(outfile)
        return merged.merged_table.view(Metadata).to_dataframe().T
    elif para == 'seq':
        merge_para = merge_seqs
    elif para == 'taxa':
        merge_para = merge_taxa
    else:
        raise(f'{para} not exists in [tab, seq, taxa]')
    if merge_para:
        merged = merge_para(qza_list)
        merged.merged_data.save(outfile)
        return merged.merged_data.view(Metadata).to_dataframe()


# %%
def merge_all_qza(jsons, outdir):
    outdict = {}
    seq_qza_list = get_qza_list(jsons, qza_name='asv_seq_qza')
    outfile = os.path.join(outdir, 'merged_seq.qza')
    merged_seq = merge_qza(seq_qza_list, outfile=outfile, para='seq')
    outdict.update({'merged_seq_qza':outfile})

    seq_qza_list = get_qza_list(jsons, qza_name='asv_tax_qza')
    outfile = os.path.join(outdir, 'merged_tax.qza')
    merged_tax = merge_qza(seq_qza_list, outfile=outfile, para='taxa')
    outdict.update({'merged_tax_qza':outfile})

    seq_qza_list = get_qza_list(jsons, qza_name='asv_tab_qza')
    outfile = os.path.join(outdir, 'merged_tab.qza')
    merged_tab = merge_qza(seq_qza_list, outfile=outfile, para='tab')
    outdict.update({'merged_tab_qza':outfile})
    
    df_table = pd.concat([merged_seq, merged_tax, merged_tab], axis=1)
    outfile = os.path.join(outdir, 'merged_tab.csv')
    df_table.to_csv(outfile)
    outdict.update({'merged_tab_csv':outfile})
    return outdict


# %%
def generate_seq_tree(seq_qza, outdir):
    from qiime2 import Artifact
    import qiime2.plugins.phylogeny.actions as phylogeny_actions
    from skbio import TreeNode
    
    rep_seqs = Artifact.load(seq_qza)
    outdict = {}
    action_results = phylogeny_actions.align_to_tree_mafft_fasttree(
        sequences=rep_seqs,
    )
    aligned_rep_seqs = action_results.alignment
    masked_aligned_rep_seqs = action_results.masked_alignment

    unrooted_tree = action_results.tree
    outfile = os.path.join(outdir, 'merged_unrooted_tree.qza')
    unrooted_tree.save(outfile)
    outdict.update({'unrooted_tree_qza':outfile})
    outfile = os.path.join(outdir, 'merged_unrooted.tree')
    unrooted_tree.view(TreeNode).write(outfile)
    outdict.update({'unrooted_tree':outfile})
    
    rooted_tree = action_results.rooted_tree
    outfile = os.path.join(outdir, 'merged_rooted_tree.qza')
    rooted_tree.save(outfile)
    outdict.update({'rooted_tree_qza':outfile})
    outfile = os.path.join(outdir, 'merged_rooted.tree')
    rooted_tree.view(TreeNode).write(outfile)
    outdict.update({'rooted_tree':outfile})
    return outdict


# %%
if __name__ == '__main__':
    bin_dir = os.path.split(os.path.realpath(__file__))[0]
    pub_path = os.path.join(bin_dir, '../pub/')
    if os.path.isdir(pub_path):
        sys.path.append(pub_path)
    else:
        raise(f'{pub_path} not exists')
    
    from write_json import write_json
    from mkdir import mkdir
    
    parse = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parse.add_argument('-i', '--input', required=True, nargs='+', help='asv dir for analysis (absolute path)')
    parse.add_argument('-o', '--outdir', required=True, help='out dir for output files')
    args = parse.parse_args()
    
    indir = args.input
    jsons = [os.path.join(i, 'Amplicon_fq2asv.json') for i in indir]
    outdir = args.outdir
    mkdir(outdir)
    logfile = os.path.join(outdir, 'log')
    logging.basicConfig(level=logging.INFO, filename=logfile, format='%(asctime)s %(levelname)s %(message)s',datefmt='%Y-%m-%d %H:%M:%S')
    
    try:
        info_dict = merge_all_qza(jsons, outdir)
    except Exception as e:
        logging.error(e)
    try:
        info_dict.update(generate_seq_tree(info_dict['merged_seq_qza'], outdir))
    except Exception as e:
        logging.error(e)
        
    json_out = write_json(info_dict, outdir=outdir)
    if not json_out:
        logging.info(f'write json failed')
        logging.info(f'{info_dict}')

# %% tags=[]
# outdirs = ['/mnt/d/Yangk/work/git/amplicon_qiime2/test_out/SRR18505774/',
#            '/mnt/d/Yangk/work/git/amplicon_qiime2/test_out/SRR18505775/',
#            '/mnt/d/Yangk/work/git/amplicon_qiime2/test_out/SRR18505770/',
#            '/mnt/d/Yangk/work/git/amplicon_qiime2/test_out/SRR18505774_silva/'
#           ]
# names = [i.split('/')[-2]  for i in outdirs]
# jsons = [os.path.join(i, 'Amplicon_fq2asv.json') for i in outdirs]

# outdir = '/mnt/d/Yangk/work/git/amplicon_qiime2/test_out/merged_out/'
# infodict = merge_all_qza(jsons, outdir)
# infodict.update(generate_seq_tree(infodict['merged_seq_qza'], outdir))

# %%

# %%

# %%
