wavemesh_executable_path=$1
dataset_path=$2
reconstructed_path=$3

${wavemesh_executable_path} test ${dataset_path} -recon ${reconstructed_path} -q 10 -l 1
python tool/decompression_ratio.py ${dataset_path} ${reconstructed_path}