#!/bin/bash
# Wraps the actual program (passed in as "$@" -- see how airfoil0.sh
# calls this script) so that each of the 4 MPI ranks on a node gets
# pinned to its own GPU, its own block of CPU memory, and its own network
# card, instead of all 4 ranks contending for shared resources.
#
# Why this matters: a MareNostrum 5 GPU node has 4 GPUs, and each GPU is
# physically closest to (faster to talk to) one specific quarter of the
# node's CPU cores/memory and one specific network interface card. If
# every MPI rank were free to use any GPU/memory/network card, they'd
# often end up using ones that are "far away" on the node's internal
# wiring, silently making every GPU transfer and network message slower.
# Pinning each rank to its own nearby set (its "NUMA domain") avoids that.
#
# OMPI_COMM_WORLD_LOCAL_RANK is set automatically by Open MPI to each
# process's rank *within this node* (0-3 here, since airfoil0.sh requests
# 4 ranks per node) -- that's what picks which case below applies to which
# process.
case ${OMPI_COMM_WORLD_LOCAL_RANK} in
0)
export CUDA_VISIBLE_DEVICES=0
export OMPI_MCA_btl_openib_if_include=mlx5_0:1
export UCX_NET_DEVICES=mlx5_0:1
numactl --membind=0 "$@"
  ;;
1)
export CUDA_VISIBLE_DEVICES=1
export UCX_NET_DEVICES=mlx5_1:1
export OMPI_MCA_btl_openib_if_include=mlx5_1:1
numactl --membind=1 "$@"
  ;;
2)
export CUDA_VISIBLE_DEVICES=2
export UCX_NET_DEVICES=mlx5_4:1
export OMPI_MCA_btl_openib_if_include=mlx5_4:1
numactl --membind=2 "$@"
  ;;
3)
export CUDA_VISIBLE_DEVICES=3
export UCX_NET_DEVICES=mlx5_5:1
export OMPI_MCA_btl_openib_if_include=mlx5_5:1
numactl --membind=3 "$@"
  ;;
esac
