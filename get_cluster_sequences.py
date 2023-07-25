from collections import defaultdict
import os
from pathlib import Path
import subprocess

import pandas as pd


def get_cluster_members_dict(mmseqs2_cluster_table: str, clusters_dir: Path) -> defaultdict:
    with open(mmseqs2_cluster_table, encoding='utf8') as cluster_table:
        if not clusters_dir.is_dir():
            clusters_dir.mkdir(parents=True, exist_ok=True)
        result = defaultdict(list)
        for line in cluster_table:
            cluster_name, cluster_member_name = line.split('\t')
            result[cluster_name].append(cluster_member_name.strip())
        return result


def write_cluster_member_ids_to_file(cluster_dict: defaultdict, clusters_dir_name: str) -> None:
    for cluster_name, cluster_members in cluster_dict.items():
        cluster_dir = Path(f'{clusters_dir_name}/{cluster_name}')
        if not cluster_dir.is_dir():
            cluster_dir.mkdir(parents=True, exist_ok=True)
        with open(f'{cluster_dir}/{cluster_name}_ids.txt', 'w', encoding='utf8') as cluster_ids:
            cluster_ids.writelines([f'{line}\n' for line in cluster_members])


def retrieve_cluster_member_sequences(clusters_dir_name: str, proteins_faa: str) -> None:
    for root, dirs, files in os.walk(clusters_dir_name):
        for directory in dirs:
            file_name = f'{root}/{directory}/{directory}_ids.txt'

            cmd = f"seqtk subseq {proteins_faa} "\
                  f"<(cut -f 1 {file_name}) "\
                  f"> {file_name.removesuffix('.txt')}.faa"
            ps = subprocess.Popen(cmd, shell=True, executable="/bin/bash",
                                  stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            _ = ps.communicate()[0]


def get_cluster_member_lengths(proteins_sizes: str) -> dict:
    with open(proteins_sizes, encoding='utf8') as proteins_sizes_file:
        result = {}
        for line in proteins_sizes_file:
            protein_id, protein_length = line.split('\t')
            result[protein_id] = protein_length.strip()
        return result


def save_cluster_member_lengths(clusters_dir_name: str, proteins_sizes_dict: dict) -> None:
    for root, dirs, files in os.walk(clusters_dir_name):
        for directory in dirs:
            with open(f'{root}/{directory}/{directory}_ids.txt', encoding='utf8') as proteins_ids_file,\
                 open(f'{root}/{directory}/{directory}_lengths.txt', 'w', encoding='utf8') as proteins_lengths_file:
                for line in proteins_ids_file:
                    protein_id = line.strip()
                    proteins_lengths_file.write(f'{protein_id}\t{proteins_sizes_dict[protein_id]}\n')


def calculate_clusters_statistics(clusters_dir_name: str, results_dir: str) -> None:
    with open(f'{results_dir}/clusters_statistics.txt', 'w', encoding='utf8') as statistics_file:
        statistics_header = 'cluster_name\tcluster_members_number\tcluster_members_mean_length'
        statistics_file.write(f'{statistics_header}\n')
        columns = ['protein_id', 'protein_length']
        for root, dirs, files in os.walk(clusters_dir_name):
            for directory in dirs:
                proteins_lengths_df = pd.read_csv(f'{root}/{directory}/{directory}_lengths.txt',
                                                  index_col=None, sep='\t',
                                                  header=None, names=columns)
                cluster_members_number = proteins_lengths_df.shape[0]
                cluster_members_mean_length = proteins_lengths_df['protein_length'].mean()
                result_line = "\t".join([directory,
                                         str(cluster_members_number),
                                         str(round(cluster_members_mean_length, 2))])

                statistics_file.write(f'{result_line}\n')


def build_mafft_alignment_for_each_cluster(clusters_dir_name: str, results_dir: str) -> None:
    for root, dirs, files in os.walk(clusters_dir_name):
        for directory in dirs:
            file_name = f'{root}/{directory}/{directory}_ids.faa'
            number_of_threads = 1

            cmd = f"mafft --thread {number_of_threads} --auto {file_name} " \
                  f"1> {file_name.removesuffix('_ids.faa')}.msa 2>> {results_dir}/mafft_stderr.txt"
            ps = subprocess.Popen(cmd, shell=True, executable="/bin/bash",
                                  stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            _ = ps.communicate()[0]


def main(mmseqs2_cluster_table: str, proteins_faa: str,
         proteins_sizes: str, results_dir: str) -> None:
    clusters_dir_name = '/'.join(mmseqs2_cluster_table.split('/')[:-2]) + '/clusters_seqs'
    clusters_dir = Path(clusters_dir_name)

    cluster_dict = get_cluster_members_dict(mmseqs2_cluster_table, clusters_dir)

    write_cluster_member_ids_to_file(cluster_dict, clusters_dir_name)

    retrieve_cluster_member_sequences(clusters_dir_name, proteins_faa)

    proteins_sizes_dict = get_cluster_member_lengths(proteins_sizes)
    save_cluster_member_lengths(clusters_dir_name, proteins_sizes_dict)

    calculate_clusters_statistics(clusters_dir_name, results_dir)

    build_mafft_alignment_for_each_cluster(clusters_dir_name, results_dir)


if __name__ == '__main__':
    mmseqs2_cluster_table = '/mnt/c/crassvirales/phylomes/results/crassvirales_refseq/' \
                            '6_protein_clustering/table_clustering.tsv'
    proteins_faa = '/mnt/c/crassvirales/phylomes/crassvirales_refseq/crassvirales.faa'
    proteins_sizes = '/mnt/c/crassvirales/phylomes/crassvirales_refseq/crassvirales_protein_sizes.txt'

    results_dir = '/mnt/c/crassvirales/phylomes/results/crassvirales_refseq/'

    main(mmseqs2_cluster_table, proteins_faa, proteins_sizes, results_dir)
