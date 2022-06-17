import sys
import os
import json
import logging

import Amplicon_merge
import Amplicon_diversity
import Amplicon_difference
import argparse

if __name__ == '__main__':
    bin_dir = os.path.split(os.path.realpath(__file__))[0]
    
    from write_json import write_json
    from mkdir import mkdir
    
    parse = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parse.add_argument('-i', '--input', required=True, nargs='+', help='asv dir for analysis (absolute path)')
    parse.add_argument('-o', '--outdir', required=True, help='out dir for output files')
    parse.add_argument('-m', '--meta', default=None, help='metadata csv ')
    parse.add_argument('-j', '--json_prefix', default='Amplicon_fq2asv', help='json file prefix')
    args = parse.parse_args()

    indir = args.input
    jsons = [os.path.join(i, f'{args.json_prefix}.json') for i in indir]
    metadata_csv = args.meta
    outdir = args.outdir
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

    if args.meta:
        try:
            Amplicon_difference.get_difference(info_dict, metadata_csv, outdir)
            info_dict.update(Amplicon_difference.get_difference_res(metadata_csv, outdir))
        except Exception as e:
            logging.error(e)
    else:
        logging.error('metadata should be given')

    json_out = write_json(info_dict, outdir=outdir)
    if not json_out:
        logging.info(f'write json failed')
        logging.info(f'{info_dict}')
