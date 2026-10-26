if [ -z "$1" ]; then
  echo "Usage: bash xcut_pipeline/run_pipeline.sh N   (1: step 1, 2: step 2, 3: steps 3-5)"
else
  MESH_NAME=oat15_slc_xcut_spline_graded_1_4_20
  OPEN=""    # set OPEN="--open" to leave out the (Lfar, 0) closure point
  C2D=./construct2d   # Construct2D executable (build it with make if missing)
  OPTIONS=oat15_slc_xcut_spline_grid_options.in   # Construct2D settings for step 2
  NELM_CLEAN=2        # step 3: degenerate elements trimmed at each closure end
  if [ -n "$OPEN" ]; then NELM_CLEAN=0; fi   # open curve: nothing to trim
  if [ "$1" == "1" ];then
    #--------STEP 1---------#
    # graded wake: growth P1 % for TE<x<2, P2 % for 2<x<3, P3 % beyond
    python3 xcut_pipeline/step1_build_curve_graded_wake.py \
        --airfoil sample_airfoils/oat15_refined.dat \
        --lfar 40.0 \
        --growth 1 1 20 \
        --x1 2.0 \
        --x2 3.0 \
        --wake-type spline \
        --pull-length 1.0 \
        $OPEN \
        --out sample_airfoils/"$MESH_NAME.dat"
    # previous step 1 (fixed 100 wake points per side, MESH_NAME=oat15_slc_xcut_spline):
    # python3 xcut_pipeline/step1_build_curve.py \
    #     --airfoil sample_airfoils/oat15_refined.dat \
    #     --lfar 40.0 --nwake 100 --wake-type spline \
    #     --out sample_airfoils/"$MESH_NAME.dat"

  elif [ "$1" == "2" ]; then
    #--------STEP 2---------#
    # Construct2D, non-interactive. It only reads settings from a file named
    # grid_options.in in the folder it runs in, so it runs in a temporary
    # folder holding a copy of $OPTIONS (any name = line removed, so the
    # output is always named after MESH_NAME) and the curve; the outputs are
    # then moved back here. Answers piped in:
    #   [y   only with --open: confirm CGRD for the open curve]
    #   GRID generate the grid, BUFF use the loaded curve as-is, QUIT
    if [ ! -x "$C2D" ]; then
      echo "Construct2D executable $C2D not found -- build it (make) or set C2D"; exit 1
    fi
    if [ ! -f "$OPTIONS" ]; then
      echo "options file $OPTIONS not found in $(pwd)"; exit 1
    fi
    C2D_ABS="$(cd "$(dirname "$C2D")" && pwd)/$(basename "$C2D")"
    RUN_DIR=$(mktemp -d)
    sed '/^[[:space:]]*[nN][aA][mM][eE][[:space:]]*=/d' "$OPTIONS" > "$RUN_DIR/grid_options.in"
    cp sample_airfoils/"$MESH_NAME.dat" "$RUN_DIR/"
    ANSWERS="GRID\nBUFF\nQUIT\n"
    if [ -n "$OPEN" ]; then ANSWERS="y\n$ANSWERS"; fi
    (cd "$RUN_DIR" && printf "$ANSWERS" | "$C2D_ABS" "$MESH_NAME.dat") > "${MESH_NAME}_construct2d.log" 2>&1
    for f in "$MESH_NAME.p3d" "${MESH_NAME}_stats.p3d" "$MESH_NAME.nmf"; do
      if [ ! -f "$RUN_DIR/$f" ]; then
        echo "Construct2D did not write $f -- see ${MESH_NAME}_construct2d.log"; rm -rf "$RUN_DIR"; exit 1
      fi
      mv "$RUN_DIR/$f" .
    done
    rm -rf "$RUN_DIR"
    grep -E "First layer wall distance|Max .*(skew|growth)|Writing" "${MESH_NAME}_construct2d.log"

  elif [ "$1" == "3" ]; then
    #--------STEP 3---------#
    python3 xcut_pipeline/step3_trim_mesh.py \
        --p3d ${MESH_NAME}.p3d --meta sample_airfoils/${MESH_NAME}.meta.json \
        --stats-p3d ${MESH_NAME}_stats.p3d \
        --nelm-clean $NELM_CLEAN \
        --out-prefix ${MESH_NAME}_trimmed
    
    #--------STEP 4---------#
    python3 xcut_pipeline/step4_build_te_box.py \
        --trimmed-p3d ${MESH_NAME}_trimmed.p3d \
        --meta ${MESH_NAME}_trimmed.meta.json \
        --out ${MESH_NAME}_te_box.p3d
    
    #--------STEP 5---------#
    python3 xcut_pipeline/step5_merge_blocks.py \
        --trimmed-p3d ${MESH_NAME}_trimmed.p3d \
        --box-p3d ${MESH_NAME}_te_box.p3d \
        --meta ${MESH_NAME}_te_box.meta.json \
        --out-prefix oat15_2block

  else
    echo "Unknown option '$1': use 1 (step 1), 2 (step 2) or 3 (steps 3-5)"; exit 1
  fi
fi  
