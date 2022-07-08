import docker
import pandas as pd
import re
import logging
import argparse
import multiprocessing
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
    image_name='wnse/qiime2:220708'
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
    image_name='wnse/qiime2:220708'
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
    image_name='wnse/qiime2:220708'
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

            try:
                shutil.copy(merge_json_file, os.path.join(outdir, f'{step}.json'))
                # os.symlink(merge_json_file, os.path.join(outdir, f'{step}.json'))
                # with open(os.path.join(outdir, f'{step}.json'),'rt') as H:
                    # info_dict.update({'ASV_info':json.load(H)})
            except Exception as e:
                logging.error(f'symlink {merge_json_file} {e}')

            s = f'{step}\tD\t'
        except Exception as e:
            logging.error(e)
            s = f'{step}\tE\t'
        try:
            write_status(status_report, s)
            post_url(taskID, step)
        except Exception as e:
            logging.error(f'{step} status {e}')


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


    # json_out = write_json(info_dict, outdir=outdir)
    # if not json_out:
        # logging.info(f'write json failed {info_dict}')
    if not args.debug:
        tmp_dir = os.path.join(outdir, 'tmp')
        try:
            shutil.rmtree(tmp_dir)
        except Exception as e:
            logging.error(e)