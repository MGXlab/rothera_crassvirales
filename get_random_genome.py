import pandas as pd


PathToFile = str


def read_ictv_vmr_table(file_name):
    df = pd.read_excel(file_name)
    df = df.fillna("unknown")
    return df


def filter_caudoviricetes(df: pd.DataFrame) -> pd.DataFrame:
    print(f'Shape of original dataframe: {df.shape}')
    caudoviricetes_df = df.query('Class == "Caudoviricetes" and Order != "Crassvirales" and Family != "unknown"')
    print(f'Shape of dataframe after filtering: {caudoviricetes_df.shape}')
    return caudoviricetes_df


def get_families_and_subfamilies_statistics(caudoviricetes_df: pd.DataFrame) -> None:
    caudoviricetes_families = caudoviricetes_df.Family.unique()
    print(f'Caudoviricetes after filtering contains: {len(caudoviricetes_families)} families')

    caudoviricetes_subfamilies = caudoviricetes_df.Subfamily.unique()
    print(f'Caudoviricetes after filtering contains: {len(caudoviricetes_subfamilies)} subfamilies')


def get_random_taxonomic_members(caudoviricetes_df: pd.DataFrame, random_number_of_genomes: int = 10) -> pd.DataFrame:
    print(f'The random number of genomes per Caudoviricetes family is equal to {random_number_of_genomes}')

    random_df = caudoviricetes_df.groupby('Family', as_index=False, group_keys=False).apply(
        lambda x: x.sample(min(random_number_of_genomes, len(x))))

    print(f'Shape of dataframe after random selection: {random_df.shape}')
    return random_df


def extract_random_taxonomic_members_from_ictv_vmr(input_file_name: PathToFile,
                                                   output_caudoviricetes_file_name: PathToFile,
                                                   output_caudoviricetes_random_file_name: PathToFile) -> None:
    df = read_ictv_vmr_table(input_file_name)
    caudoviricetes_df = filter_caudoviricetes(df)

    get_families_and_subfamilies_statistics(caudoviricetes_df)

    caudoviricetes_df.to_csv(output_caudoviricetes_file_name, sep="\t", index=False)

    random_df = get_random_taxonomic_members(caudoviricetes_df)
    random_df.to_csv(output_caudoviricetes_random_file_name, sep="\t", index=False)


if __name__ == "__main__":
    input_file_name = "/mnt/c/crassvirales/ICTV_reference_sequences/VMR_21-221122_MSL37.xlsx"
    output_caudoviricetes_file_name = "/mnt/c/crassvirales/ICTV_reference_sequences/reference_caudoviricetes_table.txt"
    output_caudoviricetes_random_file_name = "/mnt/c/crassvirales/ICTV_reference_sequences/reference_random_table.txt"

    extract_random_taxonomic_members_from_ictv_vmr(input_file_name, output_caudoviricetes_file_name,
                                                   output_caudoviricetes_random_file_name)
