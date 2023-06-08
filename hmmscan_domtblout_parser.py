#!/usr/bin/python
from typing import TextIO

from Bio import SearchIO


def start_hmmscan_domtblout_result_analysis(domtblout_file_name: str) -> None:
    result_file_name = file_name[:-4] + '_filtered_0.05.txt'
    with open(domtblout_file_name, encoding='utf8') as input_file, \
         open(result_file_name, 'w', encoding='utf8') as result_file:
        header = ('#HMM_family', 'HMM_len', 'Query_ID', 'Query_len', 'E-value', 'HMM_start', 'HMM_end',
                  'Query_start', 'Query_end', 'Coverage')
        result_file.write('\t'.join(header) + '\n')
        for qresult in SearchIO.parse(input_file, 'hmmscan3-domtab'):
            write_results_to_file(qresult, result_file)


def write_results_to_file(qresult: SearchIO._model.query.QueryResult, result_file: TextIO) -> None:
    result_line = parse_hmmscan_domtblout_result(qresult)
    if result_line is not None:
        result_file.write(f'{result_line}\n')


def parse_hmmscan_domtblout_result(qresult: SearchIO._model.query.QueryResult) -> str:
    query_id = qresult.id  # sequence ID from fasta
    query_len = qresult.seq_len
    hits = qresult.hits
    num_hits = len(hits)
    if num_hits:
        for i in range(num_hits):
            hit_evalue = hits[i].evalue  # evalue
            hmm_len = hits[i].seq_len  # target length
            hmm_aln = int(hits[i].hsps[0].hit_end) - int(hits[i].hsps[0].hit_start)  # length of alignment
            coverage = hmm_aln / float(hmm_len)  # alignment coverage
            hmm_name = hits[i].id  # target name
            if filter_hmmscan_domtblout_result(hit_evalue, coverage):
                results = (hmm_name, hmm_len, query_id, query_len, hit_evalue,
                           hits[i].hsps[0].hit_start, hits[i].hsps[0].hit_end,
                           hits[i].hsps[0].query_start, hits[i].hsps[0].query_end, coverage)
                result_line = '\t'.join(str(i) for i in results)
                return result_line


def filter_hmmscan_domtblout_result(hit_evalue: float, coverage: float,
                                    hit_evalue_threshold: float = 0.05, coverage_threshold: float = 0.3) -> bool:
    if hit_evalue < hit_evalue_threshold and coverage > coverage_threshold:
        return True
    return False


if __name__ == "__main__":
    path = '/mnt/c/crassvirales/functional_annotation/crassvirales_confirmed/hmmscan_results'
    file_names = [f'{path}/all_genomes_without_refseq_table_4_domtblout.txt',
                  f'{path}/all_genomes_without_refseq_table_11_domtblout.txt',
                  f'{path}/all_genomes_without_refseq_table_11_TGA_W_domtblout.txt',
                  f'{path}/all_genomes_without_refseq_table_11_TAG_Q_domtblout.txt',
                  f'{path}/all_genomes_without_refseq_meta_domtblout.txt',
                  f'{path}/all_genomes_without_refseq_table_15_domtblout.txt',
                  f'{path}/all_genomes_without_refseq_table_25_domtblout.txt',
                  f'{path}/all_refseq_proteins_domtblout.txt']

    for file_name in file_names:
        start_hmmscan_domtblout_result_analysis(file_name)
