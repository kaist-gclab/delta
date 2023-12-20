import sys
import os
from time import time

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print(sys.argv)
        print('usage: python evaluate_compression_ratio.py <mesh1_path> <mesh2_path>')
        sys.exit(-1)

    ply1 = sys.argv[1]
    ply2 = sys.argv[2]
    ply1_name = os.path.basename(ply1)
    ply2_name = os.path.basename(ply2)

    if not os.path.exists(ply1):
        raise FileNotFoundError(ply1)
    if not os.path.exists(ply2):
        raise FileNotFoundError(ply2)
    
    ply1_file_size = os.path.getsize(ply1)
    ply2_file_size = os.path.getsize(ply2)

    print(f'File size of {ply1_name}: {ply1_file_size} bytes')
    print(f'File size of {ply2_name}: {ply2_file_size} bytes')
    print(f'Compression ratio: {ply1_file_size / ply2_file_size:.3f}')
