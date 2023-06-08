from collections.abc import Iterable, Mapping

import pandas as pd


ProfileId = str
ProteinName = str


def get_protein_name_dict(domtblout_file_name: str,
                          profile_name_dict: Mapping[ProfileId, ProteinName]) -> dict[ProfileId, ProteinName]:
    domtblout_df = pd.read_csv(domtblout_file_name, index_col=False, sep='\t')
    hmm_family_edited_col = domtblout_df['#HMM_family'].map(profile_name_dict).fillna(domtblout_df['#HMM_family'])
    domtblout_df_edited = domtblout_df.copy()
    domtblout_df_edited['#HMM_family'] = hmm_family_edited_col

    domtblout_df_edited.to_csv(domtblout_file_name[:-4] + '_with_names.txt', sep='\t', index=False)

    domtblout_df_edited_names = domtblout_df_edited[['#HMM_family', 'Query_ID']]
    domtblout_df_edited_names_unique = domtblout_df_edited_names.drop_duplicates()
    domtblout_df_edited_names_unique.to_csv(domtblout_file_name[:-4] + '_with_names_unique.txt', sep='\t',
                                            index=False)

    protein_name_dict = dict(zip(domtblout_df_edited_names_unique['Query_ID'],
                                 domtblout_df_edited_names_unique['#HMM_family']))

    protein_name_dict = {'_'.join(k.split('_')[3:5]): v for k, v in protein_name_dict.items()}

    return protein_name_dict


def get_profile_name_dict(profile_list_file_name: str) -> dict[ProfileId, ProteinName]:
    profile_list_df = pd.read_excel(profile_list_file_name, index_col=None)
    profile_name_dict = dict(zip(profile_list_df['profile ID'], profile_list_df.nickname))
    return profile_name_dict


def make_functional_anotation_for_prodigal(annotation_file: str, annotation_file_edited: str,
                                           protein_name_dict: Mapping[ProfileId, ProteinName]) -> None:
    with open(annotation_file, encoding='utf-8') as file, \
            open(annotation_file_edited, 'w', encoding='utf-8') as result:
        for line in file:
            line = line.split('\t')
            if not line[0].startswith('#'):
                description = line[-1].split(';')
                _, protein_number = description[0].rsplit('_', 1)
                protein_name = '_'.join([line[0], protein_number])
                protein_name_function = protein_name_dict.get(protein_name, 'hp')
                protein_name_function_line = f'name={protein_name_function}'
                description_edited = ';'.join([description[0], protein_name_function_line,
                                               *description[1:]])
                line_edited = '\t'.join([*line[:-1], description_edited])

                result.write(line_edited)


def make_functional_anotation_for_refseq(annotation_file: str, annotation_file_edited: str,
                                         protein_name_dict: Mapping[ProfileId, ProteinName]) -> None:

    with open(annotation_file, encoding='utf-8') as file, \
         open(annotation_file_edited, 'w', encoding='utf-8') as result:
        for line in file:
            line = line.split('\t')
            if not line[0].startswith('#') and len(line) >= 3:
                if line[2] == 'CDS':
                    description = line[-1].split(';')
                    # 'lcl|NC_062765.1_prot_YP_010358662.1_8'
                    _, protein_name = description[0].rsplit('-', 1)
                    protein_name_function = protein_name_dict.get(protein_name, 'hp')
                    protein_name_function_line = f'name={protein_name_function}'
                    description_edited = ';'.join([description[0], protein_name_function_line,
                                                   *description[1:]])
                    line_edited = '\t'.join([*line[:-1], description_edited])

                    result.write(line_edited)


def make_functional_annotation(profile_list_file_name: str,
                               domtblout_file_names: Iterable[str], annotation_files: Iterable[str],
                               domtblout_file_name_refseq: str, annotation_file_name_refseq: str) -> None:
    profile_name_dict = get_profile_name_dict(profile_list_file_name)

    for domtblout_file_name, annotation_file in zip(domtblout_file_names, annotation_files):
        protein_name_dict = get_protein_name_dict(domtblout_file_name, profile_name_dict)
        annotation_file_edited = annotation_file[:-4] + "_edited.gff"

        make_functional_anotation_for_prodigal(annotation_file, annotation_file_edited,
                                               protein_name_dict)

    annotation_file_name_refseq_edited = annotation_file_name_refseq[:-4] + "_edited.gff"

    protein_name_dict_refseq = get_protein_name_dict(domtblout_file_name_refseq, profile_name_dict)

    make_functional_anotation_for_refseq(annotation_file_name_refseq, annotation_file_name_refseq_edited,
                                         protein_name_dict_refseq)


if __name__ == "__main__":
    profile_list_file_name = "/mnt/c/crassvirales/crassfamily_2020/profile_list.xlsx"

    hmmscan_results_path = '/mnt/c/crassvirales/functional_annotation/crassvirales_confirmed/hmmscan_results'
    domtblout_file_names = [f'{hmmscan_results_path}/all_genomes_without_refseq_meta_domtblout_filtered_0.05.txt',
                            f'{hmmscan_results_path}/all_genomes_without_refseq_table_4_domtblout_filtered_0.05.txt',
                            f'{hmmscan_results_path}/all_genomes_without_refseq_table_11_domtblout_filtered_0.05.txt',
                            f'{hmmscan_results_path}/all_genomes_without_refseq_table_11_TAG_Q_'
                            f'domtblout_filtered_0.05.txt',
                            f'{hmmscan_results_path}/all_genomes_without_refseq_table_11_TGA_W_'
                            f'domtblout_filtered_0.05.txt',
                            f'{hmmscan_results_path}/all_genomes_without_refseq_table_15_domtblout_filtered_0.05.txt',
                            f'{hmmscan_results_path}/all_genomes_without_refseq_table_25_domtblout_filtered_0.05.txt']

    domtblout_file_name_refseq = f'{hmmscan_results_path}/all_refseq_proteins_domtblout_filtered_0.05.txt'

    prodigal_annotation_path = '/mnt/c/crassvirales/functional_annotation/crassvirales_confirmed/annotations/prodigal'
    refseq_annotation_path = '/mnt/c/crassvirales/functional_annotation/crassvirales_confirmed/annotations/refseq'
    annotation_files = [f'{prodigal_annotation_path}/all_genomes_without_refseq_meta.gff',
                        f'{prodigal_annotation_path}/all_genomes_without_refseq_table_4.gff',
                        f'{prodigal_annotation_path}/all_genomes_without_refseq_table_11.gff',
                        f'{prodigal_annotation_path}/all_genomes_without_refseq_table_11_TAG_Q.gff',
                        f'{prodigal_annotation_path}/all_genomes_without_refseq_table_11_TGA_W.gff',
                        f'{prodigal_annotation_path}/all_genomes_without_refseq_table_15.gff',
                        f'{prodigal_annotation_path}/all_genomes_without_refseq_table_25.gff']

    annotation_file_name_refseq = f'{refseq_annotation_path}/all_refseq_proteins.gff'

    make_functional_annotation(profile_list_file_name,
                               domtblout_file_names, annotation_files,
                               domtblout_file_name_refseq, annotation_file_name_refseq)
