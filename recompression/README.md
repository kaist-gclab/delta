# How to build

## Install VTK 6.3

- Download VTK 6.3 and install it with CMake.

## Build this software

- Insert the line below in skeleton/CMakeLists.txt
  '''
  SET (VTK_DIR "<path_dir>")
  '''
- Note that <path_dir> is the path directory which contains VTK installed.
  - check if `VTKConfig.cmake` and `VTKConfigVersion.cmake` exist in <path_dir>
- Then build it with CMake.

## Test for 2022 SWStarLab

- After build, check 'data/xyzrgb_dragon.ply' exists and execute `test-decompression-and-recompression.sh`

## Test for 2023 SWStarLab

### How To Get Dataset

- Visit homepage '문화재청 3D 문화유산'
- Search 국보, "기마인물형", check dataset named `국보275_도기_기마인물형_뿔잔`
- Request and download the 3D model (169.4MB)
- Grayscale the model in `MeshLab` and save for test model (size decreased to 145MB)

### How to test

- After build, execute `test-decompression-and-recompression-23.sh`

```
./test-decompression-and-recompression-23.sh <executable_path> <dataset_path> <output_path_decompressed>
```
