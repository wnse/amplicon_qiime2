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
import shutil

# %%
# sys.path.append('../pub/')

# %%
import json
# import qiime2
import pandas as pd
import logging
import argparse


# %%

def import_fq(fq_list, sample_name, tmpdir='./'):
    raw_seqs = None
    try:
        from qiime2 import Artifact
        tmpmanifest = os.path.join(tmpdir, 'tmpmanifest')
        # with open(tmpmanifest, 'w') as h:
        if len(fq_list) == 2:
            df = pd.DataFrame([sample_name] + fq_list).T
            df.columns=['sample-id','forward-absolute-filepath','reverse-absolute-filepath']
            df.to_csv(tmpmanifest, index=False, sep='\t')
            raw_seqs = Artifact.import_data('SampleData[PairedEndSequencesWithQuality]', tmpmanifest, 'PairedEndFastqManifestPhred33V2')
        if len(fq_list) == 1:
            df = pd.DataFrame([sample_name]  + fq_list).T
            df.columns=['sample-id','absolute-filepath']
            df.to_csv(tmpmanifest, index=False, sep='\t')
            raw_seqs = Artifact.import_data('SampleData[SequencesWithQuality]', tmpmanifest, 'SingleEndFastqManifestPhred33V2')
    except Exception as e:
        logging.error(e)
    return raw_seqs

def get_rep_seq(raw_seqs):
    table, rep_seqs, stats = None, None, None
    try:
        import qiime2.plugins.dada2.actions as dada2_actions
        table, rep_seqs, stats = dada2_actions.denoise_single(
            demultiplexed_seqs=raw_seqs, 
            trunc_len=0
        )
    except Exception as e:
        logging.error(e)
    return table, rep_seqs, stats

def get_tax(rep_seq, ref):
    taxonomy = None
    try:
        from qiime2 import Artifact
        import qiime2.plugins.feature_classifier.actions as feature_classifier_actions
        arti_ref = Artifact.load(ref)
        taxonomy, = feature_classifier_actions.classify_sklearn(
            classifier=arti_ref,
            reads=rep_seq,
        )
    except Exception as e:
        logging.error(e)
    return taxonomy

def data2df(data):
    from qiime2 import Metadata
    return data.view(Metadata).to_dataframe()


# %% tags=[]
def stat2df(stat):
    from qiime2 import Metadata
    df_tmp = stat.view(Metadata).to_dataframe().reset_index()
    df_tmp.columns = df_tmp.columns.astype(str).str.lower().str.replace(' ','_').str.replace('-','_', regex=False).str.replace('.','_', regex=False)
    df_tmp = df_tmp.astype(str).loc[0].to_dict()
    return df_tmp

def number_describe(df_S, name=None):
    df_sta = df_S.describe()
    if name:
        df_sta.index = str(name) + df_sta.index
    return df_sta.astype(str).to_dict()

def fq2asv(fq_list, sample_name, ref, outdir):
    out_dict = {}
    try:
        raw_seqs = import_fq(fq_list, sample_name, tmpdir=outdir)
        table, rep_seqs, stats = get_rep_seq(raw_seqs)
        tax = get_tax(rep_seqs, ref)

        df_stat = stat2df(stats)
        out_dict.update({'seq_stats':df_stat})

        df_table = pd.concat([data2df(table).T.rename(columns={sample_name:'frequency'}), 
                              data2df(rep_seqs),
                              data2df(tax)],
                             axis=1)
        out_csv = os.path.join(outdir, 'asv_tax.csv')
        df_table.to_csv(out_csv)
        out_dict.update({'asv_fre_sta':number_describe(df_table['frequency'],'asv_frequency_')})
        out_dict.update({'asv_len_sta':number_describe(df_table['Sequence'].apply(len),'asv_seq_len_')})
        out_dict.update({'asv_tax_tab':out_csv})

        table_qza = os.path.join(outdir, 'asv_table.qza')
        table.save(table_qza)
        out_dict.update({'asv_table_qza':table_qza})

        rep_seqs_qza = os.path.join(outdir, 'asv_rep_seqs.qza')
        rep_seqs.save(rep_seqs_qza)
        out_dict.update({'asv_rep_seqs_qza':rep_seqs_qza})
        
        tax_qza = os.path.join(outdir, 'asv_tax.qza')
        tax.save(tax_qza)
        out_dict.update({'asv_tax_qza':tax_qza})
    except Exception as e:
        logging.error(e)
    return out_dict
#     return data2df(stats), pd.concat([data2df(table).T, data2df(rep_seqs), data2df(tax)], axis=1)


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
    parse.add_argument('-i', '--input', required=True, nargs='+', help='fastq file for analysis (absolute path)')
    parse.add_argument('-o', '--outdir', required=True, help='out dir for output files')
    parse.add_argument('-n', '--name', default='test', help='sample name')
    parse.add_argument('-r', '--ref', default=os.path.join(bin_dir, 'database/gg-13-8-99-515-806-nb-classifier.qza'))
    args = parse.parse_args()
    
    outdir = args.outdir
    mkdir(outdir)
    logfile = os.path.join(outdir, 'log')
    logging.basicConfig(level=logging.INFO, filename=logfile, format='%(asctime)s %(levelname)s %(message)s',datefmt='%Y-%m-%d %H:%M:%S')
    
    fq_list = args.input
    sample_name = args.name
    ref = args.ref
    
    info_dict = fq2asv(fq_list, sample_name, ref, outdir)
    json_out = write_json(info_dict, outdir=outdir)
    if not json_out:
        logging.info(f'write json failed')
        logging.info(f'{info_dict}')
