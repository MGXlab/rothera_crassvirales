from collections.abc import Iterable, Mapping

import pandas as pd


Color = str
ProfileId = str
ProteinName = str
PathToFile = str
FunctionalColors = Mapping[ProteinName, Color]

functional_colors = {"gene86": "#a16a2e", "PDDEXK_alpha": "#6600cc",
                     'TerL': "#8ccfb3", 'portal': "#22b2cc",
                     'gene77': "#ffb38d", 'MCP': "#0000ff",
                     'gene75': "#008000", 'gene74': "#c0d2df",
                     'gene73': "#7990b0",
                     'IHF_54': "#b07990", 'IHF_53': "#cca9b8",
                     'Ttub': "#c4c1e8", 'Tstab': "#8d89b9",
                     'gene49': "#ccbcac",
                     'primase': "#ff9900", 'SNF2': "#009900",
                     'SF1': "#74ffa2", 'DNApB': "#ac00e6",
                     'PolA': "#ff0000", 'PDDEXK_beta': "#6600cc",
                     'ATP_43b': "#009999", 'DnaB': "#ff99ff",
                     'Thy1': "#c4c1e8", 'dUTP': "#1a78ff",
                     'UDG': "#ba9b97", 'MPP': "#ffffb3",
                     'Rep_Org': "#22b2cc", 'RNR': "#cccc96",
                     'phage_O': "#008080", 'gene48b': "#aa872f",
                     'dNK': "#ffd2e5",
                     'hp': '#808080', 'other_known_functions': '#d3d3d3'}


def get_protein_name_dict(file_name: PathToFile) -> dict[ProfileId, ProteinName]:
    df = pd.read_csv(file_name, index_col=None, sep='\t')
    protein_name_dict = dict(zip(df['Query_ID'],
                                 df['#HMM_family']))

    protein_name_dict = {'_'.join(k.split('_')[3:5]): v for k, v in protein_name_dict.items()}

    return protein_name_dict


def filter_prodigal_annotation(annotation_file: PathToFile, annotation_file_edited: PathToFile,
                               protein_name_dict: Mapping[ProfileId, ProteinName],
                               functional_colors: FunctionalColors) -> None:
    with open(annotation_file, encoding='utf-8') as file, \
         open(annotation_file_edited, 'w', encoding='utf-8') as result:
        for line in file:
            line = line.split('\t')
            if not line[0].startswith('#'):
                description = line[-1].split(';')
                _, protein_number = description[0].rsplit('_', 1)
                protein_name = '_'.join([line[0], protein_number])
                protein_name_function = protein_name_dict.get(protein_name, 'hp')

                if protein_name_function not in functional_colors:
                    protein_name_function = 'other_known_functions'

                protein_name_function_line = f'name={protein_name_function}'
                description_edited = ';'.join([description[0], protein_name_function_line,
                                               *description[2:]])
                line_edited = '\t'.join([*line[:-1], description_edited])

                result.write(line_edited)


def filter_refseq_annotation(annotation_file: PathToFile, annotation_file_edited: PathToFile,
                             protein_name_dict: Mapping[ProfileId, ProteinName],
                             functional_colors: FunctionalColors) -> None:

    with open(annotation_file, encoding='utf-8') as file, \
            open(annotation_file_edited, 'w', encoding='utf-8') as result:
        for line in file:
            line = line.split('\t')
            if not line[0].startswith('#'):
                description = line[-1].split(';')
                _, protein_name = description[0].rsplit('-', 1)
                protein_name_function = protein_name_dict.get(protein_name, 'hp')

                if protein_name_function not in functional_colors:
                    protein_name_function = 'other_known_functions'

                protein_name_function_line = f'name={protein_name_function}'
                description_edited = ';'.join([description[0], protein_name_function_line,
                                               *description[2:]])
                line_edited = '\t'.join([*line[:-1], description_edited])

                result.write(line_edited)


def filter_annotations(annotation_files: Iterable[PathToFile], hmm_protein_files: Iterable[PathToFile],
                       functional_colors: FunctionalColors) -> None:
    for annotation_file, hmm_protein_file in zip(annotation_files, hmm_protein_files):
        protein_name_dict = get_protein_name_dict(hmm_protein_file)
        annotation_file_edited = annotation_file[:-4] + "_filtered.gff"
        filter_prodigal_annotation(annotation_file, annotation_file_edited,
                                   protein_name_dict, functional_colors)

    protein_name_dict = get_protein_name_dict(hmm_protein_refseq_file)
    annotation_file_edited = annotation_refseq_file[:-4] + "_filtered.gff"
    filter_refseq_annotation(annotation_refseq_file, annotation_file_edited,
                             protein_name_dict, functional_colors)


if __name__ == "__main__":
    profile_list_file_name = "/mnt/c/crassvirales/crassfamily_2020/profile_list.xlsx"

    path = '/mnt/c/crassvirales/functional_annotation/crassvirales_confirmed/annotations/prodigal'
    path_refseq = '/mnt/c/crassvirales/functional_annotation/crassvirales_confirmed/annotations/refseq'
    annotation_files = [f'{path}/all_genomes_without_refseq_meta_edited.gff',
                        f'{path}/all_genomes_without_refseq_table_4_edited.gff',
                        f'{path}/all_genomes_without_refseq_table_11_edited.gff',
                        f'{path}/all_genomes_without_refseq_table_11_TAG_Q_edited.gff',
                        f'{path}/all_genomes_without_refseq_table_11_TGA_W_edited.gff',
                        f'{path}/all_genomes_without_refseq_table_15_edited.gff',
                        f'{path}/all_genomes_without_refseq_table_25_edited.gff']

    annotation_refseq_file = f'{path_refseq}/all_refseq_proteins_edited.gff'

    path = '/mnt/c/crassvirales/functional_annotation/crassvirales_confirmed/hmmscan_results'
    hmm_protein_files = [f'{path}/all_genomes_without_refseq_meta_domtblout_filtered_0.05_with_names_unique.txt',
                         f'{path}/all_genomes_without_refseq_table_4_domtblout_filtered_0.05_with_names_unique.txt',
                         f'{path}/all_genomes_without_refseq_table_11_domtblout_filtered_0.05_with_names_unique.txt',
                         f'{path}/all_genomes_without_refseq_table_11_TAG_Q_domtblout_filtered_0.05_'
                         f'with_names_unique.txt',
                         f'{path}/all_genomes_without_refseq_table_11_TGA_W_domtblout_filtered_0.05_'
                         f'with_names_unique.txt',
                         f'{path}/all_genomes_without_refseq_table_15_domtblout_filtered_0.05_with_names_unique.txt',
                         f'{path}/all_genomes_without_refseq_table_25_domtblout_filtered_0.05_with_names_unique.txt']

    hmm_protein_refseq_file = f'{path}/all_refseq_proteins_domtblout_filtered_0.05_with_names_unique.txt'

    filter_annotations(annotation_files, hmm_protein_files, functional_colors)
