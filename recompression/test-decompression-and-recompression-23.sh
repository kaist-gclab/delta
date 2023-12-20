wavemesh_executable_path=$1
dataset_path=$2
reconstructed_path=$3
compressed_path=$4

${wavemesh_executable_path} test ${dataset_path} -recon ${reconstructed_path} -q 10 -l 1
python3 tool/evaluate_compression_ratio.py ${dataset_path} ${compressed_path}
python3 tool/evaluate_reconstruction_ratio.py ${dataset_path} ${reconstructed_path}
