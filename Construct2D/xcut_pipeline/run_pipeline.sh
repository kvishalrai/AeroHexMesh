if [ -z "$1" ]; then
  echo "Usage: Need 1 input [1->STEP1 ELSE STEP2-5]"
else
  MESH_NAME=oat15_slc_xcut_spline
  if [ "$1" == "1" ];then
    #--------STEP 1---------#
    python3 xcut_pipeline/step1_build_curve.py \
        --airfoil sample_airfoils/oat15_refined.dat \
        --lfar 80.0 --nwake 50 --wake-type spline \
        --out sample_airfoils/"$MESH_NAME.dat"

  else
    #-----------STEP 2------------#
    #----USER NEEDS TO DO THIS----#
    
    #--------STEP 3---------#
    python3 xcut_pipeline/step3_trim_mesh.py \
        --p3d ${MESH_NAME}.p3d --meta sample_airfoils/${MESH_NAME}.meta.json \
        --stats-p3d ${MESH_NAME}_stats.p3d \
        --nelm-clean 2 \
        --out-prefix ${MESH_NAME}_trimmed
    
    #--------STEP 4---------#
    python3 xcut_pipeline/step4_build_te_box.py \
        --trimmed-p3d ${MESH_NAME}_trimmed.p3d \
        --meta ${MESH_NAME}_trimmed.meta.json \
        --out ${MESH_NAME}_te_box.p3d
    
    #--------STEP 4---------#
    python3 xcut_pipeline/step5_merge_blocks.py \
        --trimmed-p3d ${MESH_NAME}_trimmed.p3d \
        --box-p3d ${MESH_NAME}_te_box.p3d \
        --meta ${MESH_NAME}_te_box.meta.json \
        --out-prefix oat15_2block
  fi
fi  
