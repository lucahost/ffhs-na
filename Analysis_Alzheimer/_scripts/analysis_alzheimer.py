import csv
import os

import networkx as nx
import pandas as pd

import measures


def filter_file_by_proteins(filename):
    filtered_data = []
    with open(filename) as proteins:
        protein_reader = csv.reader(proteins, delimiter='\t')
        for protein in protein_reader:
            if protein[0].__contains__(APP_PROTEIN) or \
                    protein[0].__contains__(HFE_PROTEIN) or \
                    protein[0].__contains__(MPO_PROTEIN) or \
                    protein[0].__contains__(NOS3_PROTEIN) or \
                    protein[0].__contains__(PLAU_PROTEIN):
                filtered_data.append(protein)

    return filtered_data


def filter_proteins_and_write_csv(source_filename, target_filename):
    print(f'start filtering file {source_filename}')
    filtered = filter_file_by_proteins(source_filename)
    print(f'found {len(filtered)} matching connections')

    write_filtered_data_to_csv(filtered, target_filename, ['Source Target neighborhood fusion cooccurence coexpression experimental database textmining combined_score'])


def write_filtered_data_to_csv(data, filename, header_row):
    with open(filename, 'w', newline='') as f:
        write = csv.writer(f)
        write.writerow(header_row)
        write.writerows(data)


def print_summary_of_separated_genes(working_dir):
    data_dir = os.path.join(working_dir, '..\\Data\\filtered_protein_data\\CSV')
    df_genes_separated = pd.concat(
        map(lambda x: pd.read_csv(x, sep=';'), [
            os.path.join(data_dir, 'APP.csv'),
            os.path.join(data_dir, 'HFE.csv'),
            os.path.join(data_dir, 'MPO.csv'),
            os.path.join(data_dir, 'NOS3.csv'),
            os.path.join(data_dir, 'PLAU.csv'),
        ]), ignore_index=True)

    print(df_genes_separated)

    G_separated = nx.from_pandas_edgelist(df_genes_separated, edge_attr='combined_score', source='Source',
                                          target='Target', create_using=nx.Graph())

    print('summary (from separated genes):')
    measures.get_summary(G_separated)
    # measures.plot_degree_histogram(G_separated)
    # nx.draw(G_separated)

def print_summary_of_combined_genes(filename):
    df_genes_separated = pd.concat(map(lambda x: pd.read_csv(x, sep='\s+'), [filename]), ignore_index=True)

    print(df_genes_separated)

    G_combined = nx.from_pandas_edgelist(df_genes_separated, edge_attr='combined_score', source='Source',
                                          target='Target', create_using=nx.Graph())

    print('summary (from separated genes):')
    measures.get_summary(G_combined)
    # measures.plot_degree_histogram(G_combined)
    nx.draw(G_combined)


if __name__ == '__main__':
    APP_PROTEIN = '9606.ENSP00000477213'
    HFE_PROTEIN = '9606.ENSP00000417404'
    MPO_PROTEIN = '9606.ENSP00000225275'
    NOS3_PROTEIN = '9606.ENSP00000297494'
    PLAU_PROTEIN = '9606.ENSP00000361850'
    PROTEINS = [APP_PROTEIN, HFE_PROTEIN, MPO_PROTEIN, NOS3_PROTEIN, PLAU_PROTEIN]

    working_dir = os.path.dirname(__file__)

    # only basic attributes
    # filename = os.path.join(working_dir, '..\\Data\\all_protein_data\\9606.protein.links.v11.5.txt')
    # target_filename = os.path.join(working_dir, '..\\Data\\all_protein_data\\9606.protein.links.v11.5_filtered.csv')
    # filter_proteins_and_write_csv(filename, target_filename)

    # with additional attributes
    filename = os.path.join(working_dir, '..\\Data\\all_proteind_data_detailed\\9606.protein.links.detailed.v11.5.txt')
    target_filename = os.path.join(working_dir, '..\\Data\\all_proteind_data_detailed\\9606.protein.links.detailed.v11.5_filtered.csv')
    # filter_proteins_and_write_csv(filename, target_filename)

    print_summary_of_separated_genes(working_dir)
    print_summary_of_combined_genes(target_filename)

