# How to build

## Install VTK 6.3

* Download VTK 6.3 and install it with CMake.  

## Build this software

* Insert the line below in skeleton/CMakeLists.txt
'''
SET (VTK_DIR "<path_dir>")
'''
* Note that <path_dir> is the path directory which contains VTK installed.
  * check if `VTKConfig.cmake` and `VTKConfigVersion.cmake` exist in <path_dir>
* Then build it with CMake.

## Test for 2022 SWStarLab

* After build, check 'data/xyzrgb_dragon.ply' exists and execute script2022.sh
