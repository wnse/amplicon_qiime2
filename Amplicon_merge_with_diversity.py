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
import logging

import Amplicon_merge
import Amplicon_diversity
import argparse

# %% tags=[]
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
    parse.add_argument('-i', '--input', required=True, nargs='+', help='asv dir for analysis (absolute path)')
    parse.add_argument('-o', '--outdir', required=True, help='out dir for output files')
    args = parse.parse_args()

    indir = args.input
    jsons = [os.path.join(i, 'Amplicon_fq2asv.json') for i in indir]
    outdir = args.outdir
    # outdir = '/mnt/d/Yangk/work/git/amplicon_qiime2/test_out/merged_out/'
    mkdir(outdir)

    logfile = os.path.join(outdir, 'log')
    logging.basicConfig(level=logging.INFO, filename=logfile, format='%(asctime)s %(levelname)s %(message)s',datefmt='%Y-%m-%d %H:%M:%S')
    info_dict = {}
    try:
        info_dict = Amplicon_merge.merge_all_qza(jsons, outdir)
    except Exception as e:
        logging.error(e)

    try:
        info_dict.update(Amplicon_merge.generate_seq_tree(info_dict['merged_seq_qza'], outdir))
    except Exception as e:
        logging.error(e)

    try:
        table_qza = info_dict['merged_tab_qza']
        tree_qza = info_dict['rooted_tree_qza']
        info_dict.update(Amplicon_diversity.get_diversity(table_qza, tree_qza, outdir))
    except Exception as e:
        logging.error(e)

    json_out = write_json(info_dict, outdir=outdir)
    if not json_out:
        logging.info(f'write json failed')
        logging.info(f'{info_dict}')
