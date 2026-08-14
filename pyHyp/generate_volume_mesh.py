"""Hyperbolic surface extrusion: a closed 3D surface mesh -> a 3D hex volume mesh.

Thin CLI wrapper around pyHyp (MDO Lab): loads pyHyp's own option set from a
JSON config file, points it at an input surface mesh, and writes the
extruded volume mesh out as Plot3D.
"""
import argparse
import json

from pyhyp import pyHyp


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, help="JSON file of pyHyp options (grid/pseudo-grid/smoothing parameters)")
    parser.add_argument("--input-file", required=True, help="Input surface mesh (CGNS or Plot3D)")
    parser.add_argument("--file-type", default="CGNS", choices=["CGNS", "PLOT3D"], help="Input surface mesh format")
    parser.add_argument("--output-file", default="volumeMesh.xyz", help="Output volume mesh (Plot3D)")
    args = parser.parse_args()

    with open(args.config) as f:
        options = json.load(f)

    options["inputFile"] = args.input_file
    options["fileType"] = args.file_type

    hyp = pyHyp(options=options)
    hyp.run()
    hyp.writePlot3D(args.output_file)


if __name__ == "__main__":
    main()
