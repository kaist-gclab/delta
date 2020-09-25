# How to build
## Install VTK 6.3
* Download VTK 6.3 and install it with CMake.  
\
## Build this software
* Insert the line below in skeleton/CMakeLists.txt
'''
SET (VTK_DIR "<path_dir>")
'''
* Note that <path_dir> is the path directory which contains VTK installed.  
* Then build it with CMake.