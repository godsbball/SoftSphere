#!/bin/bash

# Set the paths and options
NVCC=/usr/local/cuda/bin/nvcc
CCBIN=g++
INCLUDE_PATH=../../../Common
OUTPUT_OBJ=ss_beta.o
OUTPUT_BIN=ss_beta
SOURCE_FILE=ss_beta.cu

# Compilation flags
NVCC_FLAGS="-ccbin ${CCBIN} -I${INCLUDE_PATH} -m64 --threads 0 --std=c++11"
GENCODE_FLAGS="-gencode arch=compute_50,code=sm_50 \
-gencode arch=compute_52,code=sm_52 \
-gencode arch=compute_60,code=sm_60 \
-gencode arch=compute_61,code=sm_61 \
-gencode arch=compute_70,code=sm_70 \
-gencode arch=compute_75,code=sm_75 \
-gencode arch=compute_80,code=sm_80 \
-gencode arch=compute_86,code=sm_86 \
-gencode arch=compute_89,code=sm_89 \
-gencode arch=compute_90,code=sm_90 \
-gencode arch=compute_90,code=compute_90"

# Print compilation information
echo "NVCC: ${NVCC}"
echo "CCBIN: ${CCBIN}"
echo "INCLUDE_PATH: ${INCLUDE_PATH}"
echo "SOURCE_FILE: ${SOURCE_FILE}"
echo "OUTPUT_OBJ: ${OUTPUT_OBJ}"
echo "OUTPUT_BIN: ${OUTPUT_BIN}"
echo "NVCC_FLAGS: ${NVCC_FLAGS}"
echo "GENCODE_FLAGS: ${GENCODE_FLAGS}"

# Compile to object file
echo "Compiling ${SOURCE_FILE} to ${OUTPUT_OBJ}..."
${NVCC} ${NVCC_FLAGS} ${GENCODE_FLAGS} -o ${OUTPUT_OBJ} -c ${SOURCE_FILE}
echo "Compilation finished."

# Link to create the binary executable
echo "Linking ${OUTPUT_OBJ} to create ${OUTPUT_BIN}..."
${NVCC} ${NVCC_FLAGS} ${GENCODE_FLAGS} -o ${OUTPUT_BIN} ${OUTPUT_OBJ} -L/usr/local/cuda/lib64/stubs -lcuda
echo "Linking finished. Executable created: ${OUTPUT_BIN}"
