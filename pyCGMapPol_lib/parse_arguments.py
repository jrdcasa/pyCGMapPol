import argparse
import os
import datetime
import sys


# =============================================================================
def print_header(version, logger_log=None):
    msg = """
    ***********************************************************************
                        Prepare CG input files (PCIF)
              ----------------------------------------------

                                Version {}

                              Dr. Javier Ramos
                      Macromolecular Physics Department
                Instituto de Estructura de la Materia (IEM-CSIC)
                               Madrid (Spain)

        PCIF is an open-source python library to quickly  prepare the required
        files to prepare and perform CG simulations.
        Currently, only VOTCA is implemented.

        This software is distributed under the terms of the
        GNU General Public License v3.0 (GNU GPLv3). A copy of
        the license (LICENSE.txt) is included with this distribution.

    ***********************************************************************
        """.format(version)

    print(msg) if logger_log is None else logger_log.info(msg)


# =============================================================================
def command_info(opts, logger=None):

    m1 = ""
    for item in sys.argv[1:]:
        m1 += " {}".format(item)
    m = "\n\t\tCommand line: \n"
    m += "\t\t\tpython {}".format(os.path.split(sys.argv[0])[1])
    m += m1+"\n"
    m += "\t\t\t         or\n"
    m += "\t\t\tpyCGMapPol".format(os.path.split(sys.argv[0])[1])
    m += m1+"\n"
    print(m) if logger is None else logger.info(m)


# =============================================================================
def parse_arguments():
    import time

    desc = """ Prepare CG input files (PCIF).\n"""

    parser = argparse.ArgumentParser(description=desc)
    # group1 = parser.add_mutually_exclusive_group(required=True)
    parser.add_argument("-i", "--inp", dest="inpfile",
                        help="A file with information to perform the CG model.",
                        action="store", metavar="INP_FILE", required=True)

    parser.add_argument("-d" "--debug", dest="debug",
                        help="Write extra information of the AA and CG topologies.",
                        action="store_true", required=False)

    args = parser.parse_args()

    # Check for existing files:
    if not os.path.isfile(args.inpfile):
        print(desc)
        time.sleep(.25)
        parser.error("Input file must exist!!!!")

    return args
