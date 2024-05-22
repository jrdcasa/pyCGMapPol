import datetime
from pyCGMapPol_lib.CGSystemClass import CGSystemClass
from pyCGMapPol_lib.logger import init_logger
from pyCGMapPol_lib.parse_arguments import print_header, parse_arguments, command_info


def main_app(version):

    # Init logger
    logger = init_logger("Output", fileoutput="Info.log", append=False, inscreen=True)
    print_header(version, logger)

    # Parse arguments in the command line
    opts = parse_arguments()

    # Print command info
    command_info(opts, logger=logger)

    # Info
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    starttime = datetime.datetime.now()
    m = "\t\tStart Job at {} ============".format(now)
    print(m) if logger is None else logger.info(m)

    cg = CGSystemClass(log=logger, debug=opts.debug)

    m = "\t\t Read input from {}".format(opts.inpfile)
    print(m) if logger is None else logger.info(m)
    cg.read_input_from_file(opts.inpfile)

    m = "\t\t * Write XML files for mapping AA to CG models:"
    print(m) if logger is None else logger.info(m)
    cg.write_map_aa_to_cg_xml()
    m = ""
    for i in cg._xml_filename_map_aa_cg:
        m += "\t\t\t {} has been written.\n".format(i)
    print(m) if logger is None else logger.info(m)

    m = "\t\t * Write CG XML files for VOTCA.\n"
    cg.write_map_cg_xml()
    for i in cg._xml_filename_cg:
        m += "\t\t\t {} has been written.\n".format(i)
    print(m) if logger is None else logger.info(m)

    m = "\t\t * Write CG topology file for GROMACS.\n"
    cg.write_top_file()
    m += "\t\t\t topo_CG.top has been written.\n"
    print(m) if logger is None else logger.info(m)

    m = "\t\t * Write index file for GROMACS.\n"
    cg.write_index_file()
    m += "\t\t\t index_CG.ndx has been written.\n"
    print(m) if logger is None else logger.info(m)

    m = "\t\t * Write CG topology file in PSF format.\n"
    cg.write_psf_file()
    m += "\t\t\t topo_CG.psf has been written.\n"
    print(m) if logger is None else logger.info(m)

    m = "\t\t * Write map settings XML file for VOTCA.\n"
    cg.write_settings_initialdist()
    m += "\t\t\t settings_map.xml has been written.\n"
    print(m) if logger is None else logger.info(m)

    m = "\t\t * Write iterative settings XML file for VOTCA.\n"
    cg.write_settings_cg()
    m += "\t\t\t settings_cg.xml has been written.\n"
    print(m) if logger is None else logger.info(m)

    m  = "\t\tWARNING: These files serve as templates.\n"
    m2 = "\t\t         Therefore, they need to be revised and adapted to meet your requirements!!!\n"
    m1 = "\t\t"+"-"*len(m2)+"\n"
    print(m1+m+m2+m1) if logger is None else logger.info(m1+m+m2+m1)

    m  = "\t\tTRICK  : To get the coordinates for the CG map after run this program, you should.\n"
    m2 = "\t\t         executes the following command:\n"
    m2 += "\t\t         csg_map --top <TPR_FILE_FOR_AA_MODEL> --trj <TRJ_or_GRO_FILE_FOR_AA_MODEL> \n"
    m2 += "\t\t                  --cg <XML_MAP_PRODUCED> --out <GRO_FILENAME_FOR_CG> \n"
    m1 = "\t\t"+"-"*len(m)+"\n"
    print(m1+m+m2+m1) if logger is None else logger.info(m1+m+m2+m1)

    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    endtime = datetime.datetime.now()
    m = "\t\tJob  Done at {} ============\n".format(now)
    print(m) if logger is None else logger.info(m)
    m = "\t\tTotal time: {0:.2f} seconds".format((endtime-starttime).total_seconds())
    print(m) if logger is None else logger.info(m)


# =============================================================================
if __name__ == "__main__":

    __version__ = "0.1"
    main_app(__version__)
