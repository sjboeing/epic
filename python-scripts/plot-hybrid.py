#!/usr/bin/env python
import argparse
from tools.plots import plot_parcels
import os
import sys

from tools.bokeh_plots import bokeh_plot_parcels, bokeh_plot_hybrid

try:
    parser = argparse.ArgumentParser(
        description="Plot the parcels of several or individual time steps.")

    # 24 March 2021
    # https://stackoverflow.com/questions/24180527/argparse-required-arguments-listed-under-optional-arguments
    required = parser.add_argument_group('required arguments')

    required.add_argument("--filename_parcel",
                          type=str,
                          required=True,
                          help="hdf5 output file of EPIC")

    required.add_argument("--filename_field",
                          type=str,
                          required=True,
                          help="hdf5 output file of EPIC")

    parser.add_argument("--step_parcel",
                        type=int,
                        required=True,
                        help="step to plot, " \
                            + "e.g. 10")

    parser.add_argument("--step_field",
                        type=int,
                        required=True,
                        help="steps to plot, " \
                            + "e.g. 10")

    parser.add_argument("--coloring",
                        type=str,
                        required=False,
                        default='buoyancy',
                        help="how to color the ellipses")


    parser.add_argument("--fmt",
                        type=str,
                        required=False,
                        default="png",
                        help="save format (default: png)")

    parser.add_argument("--xmin",
                        type=float,
                        required=False,
                        default=None,
                        help=" float to determine x min")

    parser.add_argument("--xmax",
                        type=float,
                        required=False,
                        default=None,
                        help=" float to determine x max")

    parser.add_argument("--ymin",
                        type=float,
                        required=False,
                        default=None,
                        help="float to determine y min")

    parser.add_argument("--ymax",
                        type=float,
                        required=False,
                        default=None,
                        help="float to determine y max")

    parser.add_argument("--display",
                        type=str,
                        required=False,
                        default='full HD',
                        help="display (bokeh only)")


    if not '--filename_parcel' in sys.argv:
        parser.print_help()
        exit(0)

    if not '--filename_field' in sys.argv:
        parser.print_help()
        exit(0)

    args = parser.parse_args()

    if not os.path.exists(args.filename_parcel):
        raise IOError("File '" + args.filename_parcel + "' does not exist.")

    if not os.path.exists(args.filename_field):
        raise IOError("File '" + args.filename_field + "' does not exist.")

    bokeh_plot_hybrid(fname_parcel=args.filename_parcel,
                               fname_field=args.filename_field,
                               step_parcel=args.step_parcel,
                               step_field=args.step_field,
                               fmt=args.fmt,
                               coloring=args.coloring,
                               display=args.display,
                               xmin=args.xmin,
                               xmax=args.xmax,
                               ymin=args.ymin,
                               ymax=args.ymax)

except Exception as ex:
    print(ex)
