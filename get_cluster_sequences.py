from collections import defaultdict
import os
from pathlib import Path
import subprocess


def get_cluster_members_dict(mmseqs2_cluster_table, clusters_dir):
    with open(mmseqs2_cluster_table, encoding='utf8') as cluster_table:
        if not clusters_dir.is_dir():
            clusters_dir.mkdir(parents=True, exist_ok=True)
        result = defaultdict(list)
        for line in cluster_table:
            cluster_name, cluster_member_name = line.split('\t')
            result[cluster_name].append(cluster_member_name.strip())
        return result


def write_cluster_member_ids_to_file(cluster_dict, clusters_dir):
    for cluster_name, cluster_members in cluster_dict.items():
        cluster_dir = Path(f'{clusters_dir}/{cluster_name}')
        if not cluster_dir.is_dir():
            cluster_dir.mkdir(parents=True, exist_ok=True)
        with open(f'{cluster_dir}/{cluster_name}_ids.txt', 'w', encoding='utf8') as cluster_ids:
            cluster_ids.writelines([f'{line}\n' for line in cluster_members])


def retrieve_cluster_member_sequences(clusters_dir_name, proteins_faa):
    for root, dirs, files in os.walk(clusters_dir_name):
        for directory in dirs:
            file_name = f'{root}/{directory}/{directory}_ids.txt'

            cmd = f"seqtk subseq {proteins_faa} "\
                  f"<(cut -f 1 {file_name}) "\
                  f"> {file_name.removesuffix('.txt')}.faa"
            ps = subprocess.Popen(cmd, shell=True, executable="/bin/bash",
                                  stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            output = ps.communicate()[0]


if __name__ == '__main__':
    mmseqs2_cluster_table = '/mnt/c/crassvirales/phylomes/results/crassvirales_refseq/' \
                            '6_protein_clustering/table_clustering.tsv'
    clusters_dir_name = '/'.join(mmseqs2_cluster_table.split('/')[:-2]) + '/clusters_seqs'
    clusters_dir = Path(clusters_dir_name)
    proteins_faa = '/mnt/c/crassvirales/phylomes/crassvirales_refseq/crassvirales.faa'

    # cluster_dict = get_cluster_members_dict(mmseqs2_cluster_table, clusters_dir)

    # write_cluster_member_ids_to_file(cluster_dict, clusters_dir)

    retrieve_cluster_member_sequences(clusters_dir_name, proteins_faa)
