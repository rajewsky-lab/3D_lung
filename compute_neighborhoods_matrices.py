#!/usr/bin/env python3
__version__ = "0.1"
__author__ = [
    "Nikos Karaiskos",
]
__license__ = "GPL"
__email__ = [
    "nikolaos.karaiskos@mdc-berlin.de",
]

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
import time
from tqdm import tqdm


def compute_3d_neighborhoods_and_extract_distances():
    """
    Computes all neighborhoods and distances for every cell.
    By defauly 
    """
    # read all metadata (exported from Seurat object)
    start_time = time.time()
    print('Reading metadata table and aligned locations ...')
    md = pd.read_csv('metadata.csv', index_col=0)

    # read aligned locations for all sections and replace the unaligned in the metadata
    for section in md.section.unique():
        df = pd.read_csv(section + '_aligned.csv', header=None)

        # add index because cell_ids are missing
        df.index = md[md.section == section].index

        # change coordinates to the aligned ones
        md.loc[md.section == section, 'CenterX_global_px'] = df[0]
        md.loc[md.section == section, 'CenterY_global_px'] = df[1]

        # add z-coordinates in microns (0, 30, 60, 90, 120, 150)
        md.loc[md.section == section, 'CenterZ_global_px'] = 30*(int(section.split('_')[1])-4)//6

        # transform them to pixels
        md.loc[md.section == section, 'CenterZ_global_px'] *= 5.6


    # set neighborhood radius for 250 microns (and multiply by 5.6 for pixels)
    radius = 250 * 5.6

    # Compute neighborhoods
    neighbors = dict()

    print('Computing neighbors for', md.shape[0],'cells ...')
    for cell_i in tqdm(range(md.shape[0])):
        cell_ID_i = md.index[cell_i]

        d_i = md.iloc[[cell_i]][['CenterX_global_px', 'CenterY_global_px', 'CenterZ_global_px']]
        
        distances = cdist(d_i, md[['CenterX_global_px', 'CenterY_global_px', 'CenterZ_global_px']])
        distances = np.round(distances).astype(int)

        neighbor_cell_IDs = list(md.index[(distances <= radius).flatten()].to_numpy())
        neighbor_distances = list(distances[distances <= radius])
    
        neighbors[cell_ID_i] = [neighbor_cell_IDs, neighbor_distances]

    # Writing to CSV
    print ('Writing results to disk ...')
    with open("output/neighbors_distance_250.csv", "w") as file:
        for cell_ID, (neighbor_IDs, distances) in tqdm(neighbors.items()):
            line = ",".join([cell_ID] + neighbor_IDs + list(map(str, distances)))  # Concatenating all values
            file.write(line + "\n")  # Write line to file


def cast_distances_to_integer(line):
    """
    Helping function to parse the matrix created by the `compute_3d_neighborhoods_and_extract_distances` function.
    """
    line2 = line[:(len(line)-1) // 2 + 1] + [int(x) for x in line[((len(line)-1) // 2 + 1):]]
    return(line2)


def count_cell_types(line, md, distance):
    """
    Counts cell types found in the neighborhood for a given cell.

    Args:
        line: Parsed from the .csv output of the `compute_3d_neighborhoods_and_extract_distances` function
        md: The metadata matrix (pd.DataFrame)
        distance: The radius of the neighborhood in microns (int)
    """
    # Convert to pixels
    DISTANCE = distance * 5.6
    
    new_row = {'cell_id':line[0],
               '5_Tumor':0, '4_Tumor':0, '1_Fibroblasts':0, '3_CD8 T cells':0,
               '14_Immune proliferating':0, '2_Macrophages':0, '9_Monocytes':0,
               '13_T reg':0, '8_Endothelium vascular':0, '22_Tumor proliferating':0,
               '20_Tumor':0, '11_Pericytes':0, '7_Macrophages':0,
               '17_Endothelium lymphatic':0, '15_Basal':0, '0_Fibroblasts':0,
               '10_Plasma cells':0, '22_B Cells':0, '12_Smooth muscle':0, '23_':0,
               '19_Mast/Basophils':0, '16_Dendritic myeloid':0,
               '18_Alveolar epithelial':0, '6_Airway epithelium':0,
               '5_Tumor_min_dist':10000, '4_Tumor_min_dist':10000, '1_Fibroblasts_min_dist':10000, '3_CD8 T cells_min_dist':10000,
               '14_Immune proliferating_min_dist':10000, '2_Macrophages_min_dist':10000, '9_Monocytes_min_dist':10000,
               '13_T reg_min_dist':10000, '8_Endothelium vascular_min_dist':10000, '22_Tumor proliferating_min_dist':10000,
               '20_Tumor_min_dist':10000, '11_Pericytes_min_dist':10000, '7_Macrophages_min_dist':10000,
               '17_Endothelium lymphatic_min_dist':10000, '15_Basal_min_dist':10000, '0_Fibroblasts_min_dist':10000,
               '10_Plasma cells_min_dist':10000, '22_B Cells_min_dist':10000, '12_Smooth muscle_min_dist':10000, '23__min_dist':10000,
               '19_Mast/Basophils_min_dist':10000, '16_Dendritic myeloid_min_dist':10000,
               '18_Alveolar epithelial_min_dist':10000, '6_Airway epithelium_min_dist':10000}
    
    line = np.array(line) # cast to numpy for vectorization
    neighbor_distances = line[(((len(line)-1) // 2) + 1):]
    neighbor_distances = neighbor_distances.astype('int') # cast to int bc they're strings
    
    proximal_neighbor_idx = np.where(neighbor_distances <= DISTANCE)[0] # include only these
    
    neighbor_distances = neighbor_distances[proximal_neighbor_idx] # subset distances
    neighbor_cell_ids = line[proximal_neighbor_idx + 1] # get cell_ids
                  
    # remove own cell_id
    own_cell_id_idx = np.where(neighbor_cell_ids == line[0])[0]
    neighbor_cell_ids = np.delete(neighbor_cell_ids, own_cell_id_idx)
    neighbor_distances = np.delete(neighbor_distances, own_cell_id_idx)
        
    neighbor_annotations = md.loc[neighbor_cell_ids].annotations # get cell types
    
    annotation_counts_dict = neighbor_annotations.value_counts().to_dict() # add cell type counts 
    for cell_type in annotation_counts_dict:
        new_row[cell_type] += annotation_counts_dict[cell_type]    
    
    # find the minimum distance per cell type
    df = pd.concat([neighbor_annotations,
                    pd.DataFrame(neighbor_distances, index=neighbor_annotations.index)], axis=1)
    df_min_distance_dict = df.groupby('annotations').min().to_dict()[0]
    for cell_type in df_min_distance_dict:
        new_row[cell_type+'_min_dist'] = df_min_distance_dict[cell_type]
    
    return(new_row)    


def remove_own_cell(line):
    """
    Helping function to remove source cell from its own neighborhood.
    """
    cell_ids = pd.Series(line[:(((len(line)-1) // 2) + 1)])
    neighbor_distances = line[(((len(line)-1) // 2) + 1):]
    duplicated_idx = np.where(cell_ids.duplicated())[0]
    cell_ids = cell_ids.to_numpy()
    cell_ids = np.delete(cell_ids, duplicated_idx)
    neighbor_distances = np.delete(neighbor_distances, duplicated_idx-1)
    return(np.append(cell_ids, neighbor_distances))


def keep_only_current_section_cell_ids(line, section):
    source_cell = np.array(line[0])
    cell_ids = line[1:(((len(line)-1) // 2) + 1)]
    neighbor_distances = line[(((len(line)-1) // 2)+1):]
    
    same_section_idx = np.char.find(cell_ids.astype('U24'), section) == 0
    line = np.append(cell_ids[same_section_idx],
                     neighbor_distances[same_section_idx])
    
    return(np.append(source_cell, line))


def compute_2d_neighborhood_matrix(section, distance):
    """
    Compute 2D neighborhoods for a given section and for 50 um.

    Args:
        section: Given in the metadata format, i.e. "section_10" (str)
        distance: The radius of the neighborhood in microns (int)
    """
    print ('computing neighborhoods for', section)
    md = pd.read_csv('metadata.csv', index_col=0)
    md = md[md.section == section]

    all_dicts = []

    idx = 0
    start_time = time.time()
    with open('output/neighbors_distance_250.csv', 'r') as fi:
        for line in fi:
            idx += 1
            if idx % 5_000 == 0:
                print(f"{idx:,}", 'lines processed, total time:', round(time.time()-start_time), 'seconds')
            line = line.strip()
            line = line.split(',')
            
            source_cell = np.array(line[0])
            if np.char.find(source_cell, section) == -1:
                continue
            
            line = cast_distances_to_integer(line)
            line = remove_own_cell(line)
            line = keep_only_current_section_cell_ids(line, section)
            all_dicts.append(count_cell_types(line, md, distance))

    df = pd.DataFrame.from_records(all_dicts)
    df.to_csv('output/2D_' + section + '_50um.csv')


def compute_3d_neighborhood_matrix(distance):
    """
    Compute the full 3D neighborhood

    Args:
        distance: The radius of the neighborhood in microns (int)
    """
    md = pd.read_csv('metadata.csv', index_col=0)

    all_dicts = []

    idx = 0
    start_time = time.time()
    with open('output/neighbors_distance_250.csv', 'r') as fi:
        for line in fi:
            if idx % 5_000 == 0:
                print(f"{idx:,}", 'cells processed, total time:', round(time.time()-start_time), 'seconds')
            line = line.strip()
            line = line.split(',')
            
            line = cast_distances_to_integer(line)
            all_dicts.append(count_cell_types(line, md, distance))
            idx += 1

    df = pd.DataFrame.from_records(all_dicts)
    df.to_csv('output/3D_50um.csv')


if __name__ == "__main__":
    
    # Compute all neighborhoods and distances for 250 um
    compute_3d_neighborhoods_and_extract_distances()

    # Compute the 2D matrices for 50 um
    for section in ['section_' + str(x) for x in [10, 16, 22, 28]]:
        compute_2d_neighborhood_matrix(section, 50)

    # Compute the 3D matrix for 50 um
    compute_3d_neighborhood_matrix(50)

