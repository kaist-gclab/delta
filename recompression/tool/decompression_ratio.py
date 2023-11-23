'''
file: decompression_ratio.py
author: Yeonghun Kim
'''
import numpy as np
import pymeshlab as ml
import sys
import os
from time import time


def readPly(filename):
    ms = ml.MeshSet()
    ms.load_new_mesh(filename)
    return ms.current_mesh()


def normalize_current_mesh(ms):
    # translate, scale
    bbox = ms.current_mesh().bounding_box()
    mx, my, mz = bbox.min()
    sx, sy, sz = bbox.dim_x(), bbox.dim_y(), bbox.dim_z()
    ms.apply_filter(
        'compute_matrix_from_translation',
        axisx=-mx,
        axisy=-my,
        axisz=-mz
    )
    ms.apply_filter(
        'compute_matrix_from_scaling_or_normalization',
        uniformflag=False,
        axisx=1/sx,
        axisy=1/sy,
        axisz=1/sz,
    )


ml.pmeshlab.MeshSet.normalize_current_mesh = normalize_current_mesh

if __name__ == '__main__':

    # if len(sys.argv) < 3:
    #     print(sys.argv)
    #     print('usage: distance <mesh1_path> <mesh2_path>')
    #     sys.exit(-1)

    ply1 = sys.argv[1]
    ply2 = sys.argv[2]
    ply1_name = os.path.basename(ply1)
    ply2_name = os.path.basename(ply2)

    ms = ml.MeshSet()

    # load first mesh

    print('Loading first mesh ...')
    start = time()
    ms.load_new_mesh(ply1)
    end = time()
    print(f'{end-start:10.6f} sec elapsed. Loaded file: {ply1_name}')

    # load second mesh

    print('Loading second mesh ...')
    start = time()
    ms.load_new_mesh(ply2)
    end = time()
    print(f'{end-start:10.6f} sec elapsed. Loaded file: {ply2_name}')

    # Decompression ratio w.r.t. Hausdorff Distance

    print('Calculating decompression ratio ...')
    start = time()
    bbox_diag = ms.current_mesh().bounding_box().diagonal()
    n_vert = ms.current_mesh().vertex_number()
    value = ms.apply_filter('get_hausdorff_distance', samplenum=n_vert)
    max_hausdorff_norm = value['max'] / bbox_diag
    decompression_ratio = (1 - max_hausdorff_norm) * 100
    end = time()
    print(f'{end-start:10.6f} sec elapsed. Decompression ratio: {decompression_ratio:.3f}(%)')
