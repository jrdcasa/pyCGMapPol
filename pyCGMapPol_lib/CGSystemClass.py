import re
import os
import textwrap
from collections import defaultdict
from topology import Topology


# noinspection PyUnresolvedReferences
class CGSystemClass(object):

    __slots__ = ['_title_system', '_nkinds_aa', '_gro_filename_aa', '_top_filename_aa', '_logger',
                 '_xml_filename_map_aa_cg', '_name_topol', '_nbeads_cg', '_natoms_aa', '_nchains_aa_for_kind',
                 '_ncg_types', '_cg_dict', "_name_cg", '_debug', '_bonds_aa_dict_bykind',
                 '_mass_aa_dict_bykind', '_aggregate_index', '_type_name_bykind',
                 '_type_bead_bykind', '_bond_cg_type_dict_bykind', '_angle_cg_type_dict_bykind',
                 '_dihedral_cg_type_dict_bykind', '_atomlines_aa_dict_bykind',
                 '_cg_bond_bykind', '_bond_cg_numtables_dict', '_cg_angle_bykind',
                 '_cg_dihedral_bykind', '_angle_cg_numtables_dict', '_dihedral_cg_numtables_dict',
                 '_xml_filename_cg', '_distribution_molecules', '_map_aanamemol_to_ikind']

    # #########################################################################
    def __init__(self, debug=False, log=None):

        """
        Constructor of the class

        :param log:
        """

        self._title_system = None
        self._nkinds_aa = None
        self._gro_filename_aa = None
        self._top_filename_aa = None
        self._debug = debug

        self._name_cg = []
        self._xml_filename_map_aa_cg = []
        self._xml_filename_cg = []
        self._name_topol = []

        self._nbeads_cg = defaultdict()
        self._natoms_aa = defaultdict()
        self._nchains_aa_for_kind = defaultdict()
        self._ncg_types = defaultdict()
        self._aggregate_index = defaultdict(list)

        self._cg_dict = defaultdict(dict)

        self._bonds_aa_dict_bykind = defaultdict()
        self._mass_aa_dict_bykind = defaultdict()
        self._type_bead_bykind = defaultdict()
        self._type_name_bykind = defaultdict()
        self._cg_bond_bykind = defaultdict()
        self._bond_cg_type_dict_bykind = defaultdict()
        self._cg_angle_bykind = defaultdict()
        self._angle_cg_type_dict_bykind = defaultdict()
        self._cg_dihedral_bykind = defaultdict()
        self._dihedral_cg_type_dict_bykind = defaultdict()
        self._atomlines_aa_dict_bykind = defaultdict()
        self._bond_cg_numtables_dict = defaultdict()
        self._angle_cg_numtables_dict = defaultdict()
        self._dihedral_cg_numtables_dict = defaultdict()
        self._distribution_molecules = list()
        self._map_aanamemol_to_ikind = defaultdict()

        self._logger = log

    # #########################################################################
    # noinspection PyUnresolvedReferences
    def read_input_from_file(self, finpname):

        """
        Read the input file where the required parameters
        are stored in order to run the model CG.

        :param finpname:
        :return:
        """

        with open(finpname, 'r') as finp:
            # Read all lines
            contents = finp.readlines()
            # Remove lines with comments and return carriage
            contents = [iline.strip("\n") for iline in contents
                        if not (re.match(r"^(\s*)#", iline) or len(iline.strip()) == 0)]

            # ================= GENERAL ========================
            try:
                start_idx = contents.index("<GENERAL>")
                end_idx = contents.index("</GENERAL>")
            except ValueError as e:
                m = "\n\t\tERROR: {}\n".format(e)
                m += "\t\t\t <GENERAL> ... </GENERAL> labels must exist in the file:\n"
                m += "{}".format(finpname)
                print(m) if self._logger is None else self._logger.error(m)
                exit()

            self._title_system = contents[start_idx+1]
            self._nkinds_aa = contents[start_idx+2]
            self._gro_filename_aa = contents[start_idx+3].strip()
            self._top_filename_aa = contents[start_idx+4].strip()

            self._check_parameters("general")

            # ================= NAMES ========================
            try:
                start_idx = contents.index("<NAMES>")
                end_idx = contents.index("</NAMES>")
            except ValueError as e:
                m = "\n\t\tERROR: {}\n".format(e)
                m += "\t\t\t <NAMES> ... </NAMES> labels must exist in the file:\n"
                m += "{}".format(finpname)
                print(m) if self._logger is None else self._logger.error(m)
                exit()

            for idx in range(start_idx+1, end_idx):
                iline = contents[idx].split()
                self._name_topol.append(iline[0])
                self._name_cg.append(iline[1])
                self._xml_filename_map_aa_cg.append(iline[1] + ".xml")
                self._xml_filename_cg.append(iline[1] + "_CG.xml")

            self._check_parameters("names")

            # ================= MAP ========================
            try:
                start_idx = contents.index("<MAPS>")
                end_idx = contents.index("</MAPS>")
            except ValueError as e:
                m = "\n\t\tERROR: {}\n".format(e)
                m += "\t\t\t <MAPS> ... </MAPS> labels must exist in the file:\n"
                m += "{}".format(finpname)
                print(m) if self._logger is None else self._logger.error(m)
                exit()

            ikind = 0
            for idx in range(start_idx+1, end_idx):
                iline = contents[idx].split()
                self._natoms_aa[ikind] = iline[0]
                self._nbeads_cg[ikind] = iline[1]
                self._nchains_aa_for_kind[ikind] = iline[2]
                ikind += 1

            self._check_parameters("maps")

            # Read information from the atomistic model
            self._get_aa_topology()

            # ================= BEADS ========================
            try:
                start_idx = contents.index("<BEADS>")
                end_idx = contents.index("</BEADS>")
            except ValueError as e:
                m = "\n\t\tERROR: {}\n".format(e)
                m += "\t\t\t <BEADS> ... </BEADS> labels must exist in the file:\n"
                m += "{}".format(finpname)
                print(m) if self._logger is None else self._logger.error(m)
                exit()

            self._check_parameters("beads_previous", contents, start_idx, end_idx)

            idx = start_idx + 1
            while idx < end_idx:
                ikind = int(contents[idx].split()[0])
                type_cg_dict = defaultdict(dict)
                self._ncg_types[ikind] = int(contents[idx].split()[1])
                idx += 1
                for _ in range(self._ncg_types[ikind]):
                    iline = contents[idx].split()
                    type_cg = iline[0]
                    type_cg_dict[type_cg] = {'natoms': None, 'weights': [],
                                             'ini_idx_bead': None, 'end_idx_bead': None,
                                             'idx_aa_list': []}
                    natoms_aa_per_cg = int(iline[1])
                    weight_indexes = [float(i) for i in iline[2:-1]]
                    ktoken = iline[-1]
                    type_cg_dict[type_cg]['natoms'] = natoms_aa_per_cg
                    type_cg_dict[type_cg]['weights'] = weight_indexes
                    idx += 1

                    # Find atomistic indexes for each CG beads
                    # Use the topology_aa_dict graphs
                    ini_idx_bead = int(ktoken.split("-")[0])
                    try:
                        end_idx_bead = int(ktoken.split("-")[1])
                    except IndexError:
                        end_idx_bead = ini_idx_bead
                    type_cg_dict[type_cg]['ini_idx_bead'] = ini_idx_bead
                    type_cg_dict[type_cg]['end_idx_bead'] = end_idx_bead

                self._cg_dict[ikind] = type_cg_dict

            self._check_parameters("beads")

            # ================= INDEXES ========================
            try:
                start_idx = contents.index("<INDEXES>")
                end_idx = contents.index("</INDEXES>")
            except ValueError as e:
                m = "\n\t\tERROR: {}\n".format(e)
                m += "\t\t\t <INDEXES> ... </INDEXES> labels must exist in the file:\n"
                m += "{}".format(finpname)
                print(m) if self._logger is None else self._logger.error(m)
                exit()

            idx = start_idx + 1
            while idx < end_idx:
                type_name_tmp = []
                type_bead_tmp = []
                ikind = int(contents[idx].split()[0])
                token2 = int(contents[idx].split()[1])
                idx += 1
                idx_local = 0
                if token2 != self._ncg_types[ikind]:
                    m = "\n\t\tERROR: {}\n".format(e)
                    m += "\t\t\t <INDEXES> ... </INDEXES> is not correct:\n"
                    m += "{}".format(finpname)
                    print(m) if self._logger is None else self._logger.error(m)
                    exit()
                for _ in range(self._ncg_types[ikind]):
                    token1 = contents[idx].split()[0]
                    indices = [int(i) for i in contents[idx].split()[1:]]
                    nat_per_bead = self._cg_dict[ikind][token1]['natoms']
                    ibead = self._cg_dict[ikind][token1]['ini_idx_bead']
                    jbead = self._cg_dict[ikind][token1]['end_idx_bead']
                    weights = self._cg_dict[ikind][token1]['weights']
                    for idx_local_bead in range(0, jbead-ibead+1):
                        if any(i < 1.0 for i in weights):
                            idx_tmp_list = [i + idx_local_bead * (nat_per_bead-1) for i in indices]
                        else:
                            idx_tmp_list = [i + idx_local_bead * nat_per_bead for i in indices]
                        type_name_tmp.append(token1)
                        type_bead_tmp.append("{}{}".format(token1, idx_local+1))
                        self._cg_dict[ikind][token1]['idx_aa_list'].append(idx_tmp_list)
                        self._aggregate_index[ikind].append(idx_tmp_list)
                        idx_local += 1
                    idx += 1
                self._type_bead_bykind[ikind] = type_name_tmp
                self._type_name_bykind[ikind] = type_bead_tmp

        # Create CG topology from the information of the atomistic model
        self._create_cg_topology()

    # #########################################################################
    def write_map_aa_to_cg_xml(self):

        """
        Write the map for each CG model found in the AA topology.

        :return:
        """

        # # Read AA-GRO coordinates ==================================
        # with open(self._gro_filename_aa, 'r') as fgro:
        #     atom_contents = fgro.readlines()[2:-1]

        # WRITE XML with mapping
        # NAME OF CG REPRESENTATION ================================
        for ikind in range(self._nkinds_aa):

            with open(self._xml_filename_map_aa_cg[ikind], 'w') as fout:

                fout.writelines("<cg_molecule>\n")

                line = "  <name>{0:s}</name>\n".format(self._name_cg[ikind])
                fout.writelines(line)

                # NAME OF THE MOLECULE IN THE TOP FILE (AA REPRESENTATION)
                line = "  <ident>{0:s}</ident>\n".format(self._name_topol[ikind])
                fout.writelines(line)

                # TOPOLOGY --> DEFINING THE CG BEADS
                fout.writelines("  <topology>\n")
                fout.writelines("   <cg_beads>\n")

                # WRITE CG BEADS
                # idx_global = 0
                for ibead in range(0, self._nbeads_cg[ikind]):
                    fout.writelines("     <cg_bead>\n")
                    fout.writelines("       <name>{}</name>\n".format(self._type_name_bykind[ikind][ibead]))
                    fout.writelines("       <type>{}</type>\n".format(self._type_bead_bykind[ikind][ibead]))
                    fout.writelines("       <mapping>{}</mapping>\n".format(self._type_bead_bykind[ikind][ibead]))

                    string_beads = ""
                    for idx in self._aggregate_index[ikind][ibead]:
                        # iline = atom_contents[idx_global]
                        iline = self._atomlines_aa_dict_bykind[ikind][idx].split()
                        resid = iline[2]
                        resname = iline[3]
                        atomname = iline[4]
                        string_beads += resid.strip() + ":" + resname.strip() + ":" + atomname.strip() + " "

                    fout.writelines("       <beads>{}</beads>\n".format(string_beads))
                    fout.writelines("     </cg_bead>\n")

                fout.writelines("   </cg_beads>\n")

                # BONDING OF THE CG MODEL
                fout.writelines("   <cg_bonded>\n")
                # BONDS
                for type_bond, bond_names in self._bond_cg_type_dict_bykind[ikind].items():
                    fout.writelines("     <bond>\n")
                    fout.writelines("       <name>{}</name>\n".format(type_bond))
                    fout.writelines("     <beads>\n")
                    for ibond in bond_names:
                        fout.writelines("       {}\n".format(ibond))
                    fout.writelines("     </beads>\n")
                    fout.writelines("     </bond>\n")

                # ANGLES
                for type_angle, angle_names in self._angle_cg_type_dict_bykind[ikind].items():
                    fout.writelines("     <angle>\n")
                    fout.writelines("       <name>{}</name>\n".format(type_angle))
                    fout.writelines("     <beads>\n")
                    for iangle in angle_names:
                        fout.writelines("       {}\n".format(iangle))
                    fout.writelines("     </beads>\n")
                    fout.writelines("     </angle>\n")

                # DIHEDRALS
                for type_dihedral, dihedral_names in self._dihedral_cg_type_dict_bykind[ikind].items():
                    fout.writelines("     <dihedral>\n")
                    fout.writelines("       <name>{}</name>\n".format(type_dihedral))
                    fout.writelines("     <beads>\n")
                    for idih in dihedral_names:
                        fout.writelines("       {}\n".format(idih))
                    fout.writelines("     </beads>\n")
                    fout.writelines("     </dihedral>\n")

                fout.writelines("   </cg_bonded>\n")
                fout.writelines("  </topology>\n")

                # MAPS
                fout.writelines("  <maps>\n")
                for ibead, ivalues in self._cg_dict[ikind].items():
                    fout.writelines("    <map>\n")
                    fout.writelines("      <name>{}</name>\n".format(ibead))
                    fout.writelines("      <weights>")
                    local = 0
                    mtotal = 0
                    for idx in ivalues['idx_aa_list'][0]:
                        fout.writelines("{} ".format(self._mass_aa_dict_bykind[ikind][idx]))
                        mtotal += self._mass_aa_dict_bykind[ikind][idx]
                        local += 1
                    self._cg_dict[ikind][ibead]["mass_bead"] = mtotal
                    fout.writelines("</weights>\n".format(ibead))
                    fout.writelines("    </map>\n")
                fout.writelines("  </maps>\n")

                fout.writelines("</cg_molecule>\n")

    # #########################################################################
    def write_top_file(self, title="Polymer", nrexcl=3):

        """
        Write GROMACS topology CG file

        :param title:
        :param nrexcl:
        :return:
        """

        with open("topo_CG.top", 'w') as ftop:

            # Title section
            ftop.writelines(";\n")
            ftop.writelines(";                        CG Topology file\n")
            ftop.writelines(";                        ----------------\n")
            ftop.writelines("\n")
            ftop.writelines(";                     {} CG Map 1\n".format(title))
            ftop.writelines(";\n")

            # [defaults] section
            ftop.writelines("[ defaults ]\n")
            ftop.writelines("; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ\n")
            ftop.writelines("    1             1               no               0.0    0.0\n")
            ftop.writelines("\n")

            # [atomtypes] section
            ftop.writelines("[ atomtypes ]\n")
            ftop.writelines("; type     mass   charge   ptype  sigma  epsilon\n")
            for ikind in range(self._nkinds_aa):
                for ikey, ivalues in self._cg_dict[ikind].items():
                    ftop.writelines("{0:>4s}   {1:9.5}  0.0000  A  1.00   1.00\n".format(ikey, ivalues["mass_bead"]))
            ftop.writelines("\n")

            for ikind in range(self._nkinds_aa):

                # [ moleculetype ] section
                # pol_name = "sPPCG"
                # nrexcl = 3
                ftop.writelines("[ moleculetype ]\n")
                ftop.writelines("; Name            nrexcl\n")
                ftop.writelines("%s  %d\n" % (self._name_cg[ikind], nrexcl))
                ftop.writelines("\n")

                # [ atoms ] section
                ftop.writelines("[ atoms ]\n")
                ftop.writelines(";   nr       type  resnr residue  atom   cgnr     charge       "
                                "mass  typeB    chargeB      massB\n")
                index = 1
                for jbead in self._type_bead_bykind[ikind]:
                    ires = 1
                    charge = 0.0
                    itype = jbead
                    nameres = self._name_cg[ikind]
                    iname = self._type_name_bykind[ikind][index - 1]
                    mass = self._cg_dict[ikind][jbead]["mass_bead"]
                    iiline = ("{0:<8d} {1:<4s} {2:<5d} {3:<3s}   {4:<6s}  {5:>8d} {6:<5.3f} {7:<f}\n".format
                              (index, itype, ires, nameres, iname, index, charge, mass))
                    ftop.writelines(iiline)
                    index += 1
                ftop.writelines("\n")

                # [ bonds ] section
                ftop.writelines("[ bonds ]\n")
                ftop.writelines(";  ai    aj funct           table   k\n")
                for ibond in self._cg_bond_bykind[ikind]:
                    iat = self._type_bead_bykind[ikind][ibond[0]]
                    jat = self._type_bead_bykind[ikind][ibond[1]]
                    ilabel_bond = "-".join(sorted([iat, jat]))
                    itable = self._bond_cg_numtables_dict[ilabel_bond]
                    ik = 1.0
                    iline = ("{0:<6d} {1:<6d} {2:<6d} {3:<6d} {4:<6.1f}\n".
                             format(ibond[0] + 1, ibond[1] + 1, 8, itable, ik))
                    ftop.writelines(iline)
                ftop.writelines("\n")

                # [ angles ] section
                ftop.writelines("[ angles ]\n")
                ftop.writelines(";  ai    aj  ak table   k          \n")
                for iangle in self._cg_angle_bykind[ikind]:
                    iat = self._type_bead_bykind[ikind][iangle[0]]
                    jat = self._type_bead_bykind[ikind][iangle[1]]
                    kat = self._type_bead_bykind[ikind][iangle[2]]
                    ill = sorted([iat, kat])
                    ilabel_angle = ill[0] + "-" + jat + "-" + ill[1]
                    itable = self._angle_cg_numtables_dict[ilabel_angle]
                    ik = 1.0
                    iline = ("{0:<6d} {1:<6d} {2:<6d} {3:<6d} {4:<6d} {5:<6.1f}\n"
                             .format(iangle[0] + 1, iangle[1] + 1, iangle[2] + 1, 8, itable, ik))
                    ftop.writelines(iline)
                ftop.writelines("\n")

                # [ dihedrals ] section
                ftop.writelines("[ dihedrals ]\n")
                ftop.writelines(";  ai    aj  ak   al table   k          \n")
                for itors in self._cg_dihedral_bykind[ikind]:
                    iat = self._type_bead_bykind[ikind][itors[0]]
                    jat = self._type_bead_bykind[ikind][itors[1]]
                    kat = self._type_bead_bykind[ikind][itors[2]]
                    lat = self._type_bead_bykind[ikind][itors[3]]
                    ilabel_dihedral = iat + "-" + jat + "-" + kat + "-" + lat
                    ilabel_dihedral_inv = lat + "-" + kat + "-" + jat + "-" + iat
                    if ilabel_dihedral in self._dihedral_cg_numtables_dict.keys():
                        itable = self._dihedral_cg_numtables_dict[ilabel_dihedral]
                    elif ilabel_dihedral_inv in self._dihedral_cg_numtables_dict.keys():
                        itable = self._dihedral_cg_numtables_dict[ilabel_dihedral_inv]
                    else:
                        print("ERROR: Dihedral table not found. Label {} or {}".
                              format(ilabel_dihedral, ilabel_dihedral_inv))
                        exit()
                    ik = 1.0
                    iline = ("{0:<6d} {1:<6d} {2:<6d} {3:<6d} {4:<6d} {5:<6d} {6:<6.1f}\n".
                             format(itors[0] + 1, itors[1] + 1, itors[2] + 1,
                                    itors[3] + 1, 8, itable, ik))
                    ftop.writelines(iline)
                ftop.writelines("\n")

            # [ System ] section
            ftop.writelines("[ system ]\n")
            ftop.writelines("{}\n".format(self._title_system))
            ftop.writelines("\n")

            # [ Molecules ] section
            ftop.writelines("[ molecules ]\n")
            ftop.writelines("; Compound        #mols\n")
            for ikind in range(self._nkinds_aa):
                ftop.writelines("%s  %d\n" % (self._name_cg[ikind], self._nchains_aa_for_kind[ikind]))
            ftop.writelines("\n")

    # #########################################################################
    def write_index_file(self):

        """
        Write an index file in a format ndx for GROMACS
        :return:
        """

        d = defaultdict(list)
        index = 1
        with open("index_CG.ndx", 'w') as fout:
            imon = 0
            fout.writelines("[ System ]\n")
            write_lines = ""
            for ikind in range(self._nkinds_aa):
                nbeads = len(self._type_bead_bykind[ikind])
                for ich in range(self._nchains_aa_for_kind[ikind]):
                    for _ in range(nbeads):
                        n = (imon+1)
                        write_lines += "{0:<6d} ".format(n)
                        imon += 1
                for ich in range(self._nchains_aa_for_kind[ikind]):
                    for b_name in self._type_bead_bykind[ikind]:
                        d[b_name].append(index)
                        index += 1

            wrapper = textwrap.TextWrapper(width=60)
            word_list = wrapper.wrap(text=write_lines)
            for element in word_list:
                fout.writelines(element)
                fout.writelines("\n")


            for i in d:
                write_lines = ""
                fout.writelines("[ %s ]" % i)
                for j in d[i]:
                    write_lines += "{0:<6d} ".format(int(j))
                fout.writelines("\n")
                wrapper = textwrap.TextWrapper(width=60)
                word_list = wrapper.wrap(text=write_lines)
                for element in word_list:
                    fout.writelines(element)
                    fout.writelines("\n")


            ich = 0
            imon = 0
            for ikind in range(self._nkinds_aa):
                for _ in range(self._nchains_aa_for_kind[ikind]):
                    write_lines = ""
                    fout.writelines("[ mol{0:d} ]".format(ich+1))
                    ich += 1
                    for _ in range(self._nbeads_cg[ikind]):
                        n = imon + 1
                        write_lines += "{0:<6d} ".format(int(n))
                        imon += 1
                    fout.writelines("\n")
                    wrapper = textwrap.TextWrapper(width=60)
                    word_list = wrapper.wrap(text=write_lines)
                    for element in word_list:
                        fout.writelines(element)
                        fout.writelines("\n")

    # #########################################################################
    def write_psf_file(self):
        """
        Write CG topology in PSF format

        :return:
        """

        with open("topo_CG.psf", 'w') as f:

            # Title section
            f.writelines("PSF\n")
            f.writelines("\n")
            iiline = str("%8d %s\n" % (1, '!NTITLE'))
            f.writelines(iiline)
            # TODO
            iiline = " !TODO: title PSF\n"
            f.writelines(iiline)
            f.writelines("\n")

            # Natoms section
            ncg_beads = 0
            nbonds = 0
            ntheta = 0
            nphi = 0
            for ikind in range(self._nkinds_aa):
                ncg_beads += self._nchains_aa_for_kind[ikind] * self._nbeads_cg[ikind]
                nbonds += self._nchains_aa_for_kind[ikind] * len(self._cg_bond_bykind[ikind])
                ntheta += self._nchains_aa_for_kind[ikind] * len(self._cg_angle_bykind[ikind])
                nphi += self._nchains_aa_for_kind[ikind] * len(self._cg_dihedral_bykind[ikind])

            iiline = str("{0:8d} {1:s}\n".format(ncg_beads, '!NATOM'))
            f.writelines(iiline)

            iglobal = 0
            ich = 0
            for item in self._distribution_molecules:
                ikind = self._map_aanamemol_to_ikind[item[0]]
                inmols = item[1]
                for _ in range(inmols):
                    ilocal = 0
                    for _ in self._type_bead_bykind[ikind]:
                        typechain = str(ikind+1)
                        tmp = 0
                        charge = 0.0
                        mass = 14.0
                        iatname = self._type_name_bykind[ikind][ilocal]
                        iattype = self._type_bead_bykind[ikind][ilocal]
                        nameres = self._name_cg[ikind][0:4]
                        iiline = str("{0:8d} {1:<4s} {2:<4d} {3:>4s} {4:<4s} {5:<4s} {6:<14.6f}{7:<14.6f}{8:8d}\n"
                                     .format(iglobal + 1, typechain, ich+1, nameres, iatname[0:4],
                                             iattype[0:4], charge, mass, tmp))
                        iglobal += 1
                        ilocal += 1
                        f.writelines(iiline)
                    ich += 1
            f.writelines("\n")
            f.flush()

            # NBonds section ======================================================
            iline = str("{0:8d} {1:s}\n".format(nbonds, '!NBOND: bonds'))
            f.writelines(iline)
            idx_line = 1
            acc = 0
            for item in self._distribution_molecules:
                ikind = self._map_aanamemol_to_ikind[item[0]]
                inmols = item[1]
                natch = len(self._type_name_bykind[ikind])
                ich = 0
                for _ in range(inmols):
                        for ib in self._cg_bond_bykind[ikind]:
                            iat0 = ib[0]+1+ich*natch + acc
                            iat1 = ib[1]+1+ich*natch + acc
                            iiline = str("{0:8d}{1:8d}").format(iat1, iat0)
                            if idx_line % 4 == 0:
                                f.writelines(iiline)
                                f.writelines("\n")
                            else:
                                f.writelines(iiline)
                            idx_line += 1
                        ich += 1
                acc += inmols * natch
            if nbonds % 4 != 0:
                f.writelines("\n\n")
            else:
                f.writelines("\n")

            # NAngles section =====================================================
            line0 = str("{0:8d} {1:s}\n".format(ntheta, '!NTHETA: angles'))
            f.writelines(line0)
            idx_line = 1
            acc = 0
            for item in self._distribution_molecules:
                ikind = self._map_aanamemol_to_ikind[item[0]]
                inmols = item[1]
                natch = len(self._type_name_bykind[ikind])
                ich = 0
                for _ in range(inmols):
                    for ia in self._cg_angle_bykind[ikind]:
                        iat0 = ia[0]+1+ich*natch + acc
                        iat1 = ia[1]+1+ich*natch + acc
                        iat2 = ia[2]+1+ich*natch + acc
                        iiline = str("{0:8d}{1:8d}{2:8d}").format(iat2, iat1, iat0)
                        if idx_line % 3 == 0:
                            f.writelines(iiline)
                            f.writelines("\n")
                        else:
                            f.writelines(iiline)
                        idx_line += 1
                    ich += 1
                acc += inmols * natch
            if ntheta % 3 != 0:
                f.writelines("\n\n")
            else:
                f.writelines("\n")

            # NPhi section ========================================================
            line0 = str("{0:8d} {1:s}\n".format(nphi, '!NPHI: dihedrals'))
            f.writelines(line0)
            acc = 0
            idx_line = 1
            for item in self._distribution_molecules:
                ikind = self._map_aanamemol_to_ikind[item[0]]
                inmols = item[1]
                natch = len(self._type_name_bykind[ikind])
                ich = 0
                for _ in range(inmols):
                    for iphi in self._cg_dihedral_bykind[ikind]:
                        iat0 = iphi[0]+1+ich*natch + acc
                        iat1 = iphi[1]+1+ich*natch + acc
                        iat2 = iphi[2]+1+ich*natch + acc
                        iat3 = iphi[3]+1+ich*natch + acc
                        iiline = str("{0:8d}{1:8d}{2:8d}{3:8d}").format(iat3, iat2, iat1, iat0)
                        if idx_line % 2 == 0:
                            f.writelines(iiline)
                            f.writelines("\n")
                        else:
                            f.writelines(iiline)
                        idx_line += 1
                    ich += 1
                acc += inmols * natch
            if nphi % 2 != 0:
                f.writelines("\n\n")
            else:
                f.writelines("\n")

            # NImp section ========================================================
            line0 = str("{0:8d} {1:s}\n".format(0, '!NIMPHI: impropers'))
            f.writelines(line0)
            f.writelines("\n\n\n")

            line0 = str("{0:8d} {1:s}\n".format(0, '!NDON: donors'))
            f.writelines(line0)
            f.writelines("\n\n\n")
            line0 = str("{0:8d} {1:s}\n".format(0, '!NACC: acceptors'))
            f.writelines(line0)
            f.writelines("\n\n\n")
            line0 = str("{0:8d} {1:s}\n".format(0, '!NNB'))
            f.writelines(line0)
            f.writelines("\n\n\n")
            line0 = str("{0:8d}{1:8d} {2:s}\n".format(1, 0, '!NGRP: acceptors'))
            f.writelines(line0)
            f.writelines("\n\n\n")
            line1 = str("{0:8d}{1:8d}{2:8d}\n".format(0, 0, 0))
            f.writelines(line1)

            return

    # #########################################################################
    def _get_aa_topology(self):

        # Deduce CG bonds from AA bonds ============================
        with open(self._top_filename_aa, 'r') as ftop:
            topo_contents = ftop.readlines()
            indices_lines = defaultdict(list)
            for num, line in enumerate(topo_contents, 1):
                if re.match(r"^\[(\s*)bonds(\s*)\]", line):
                    indices_lines["bonds"].append(num - 1)
                if re.match(r"^\[(\s*)atoms(\s*)\]", line):
                    indices_lines["atoms"].append(num - 1)
                if re.match(r"^\[(\s*)molecules(\s*)\]", line):
                    indices_lines["molecules"].append(num - 1)
                if re.match(r"^\[(\s*)moleculetype(\s*)\]", line):
                    indices_lines["moleculetype"].append(num - 1)

            for ikind in range(0, self._nkinds_aa):
                bonds_aa_list = []
                mass_aa_list = []
                atomlines_aa_dict = []
                bonds_aa_dict = defaultdict(list)
                idx = indices_lines["bonds"][ikind] + 1
                while True:
                    iline = topo_contents[idx]
                    if re.match(r"^(\s*);", iline):
                        idx += 1
                        continue
                    try:
                        i, j, typebond, k, b, = iline.split()
                        ibond = int(i) - 1
                        jbond = int(j) - 1
                        bonds_aa_list.append([ibond, jbond])
                        bonds_aa_dict[ibond].append(int(jbond))
                        bonds_aa_dict[jbond].append(ibond)
                        idx += 1
                    except ValueError:
                        break

                idx = indices_lines["atoms"][ikind] + 1
                while True:
                    iline = topo_contents[idx]
                    if re.match(r"^(\s*);", iline):
                        idx += 1
                        continue
                    try:
                        tokens = iline.split()
                        mass_aa_list.append(float(tokens[7]))
                        atomlines_aa_dict.append(iline)
                        idx += 1
                    except IndexError:
                        break
                self._atomlines_aa_dict_bykind[ikind] = atomlines_aa_dict
                self._bonds_aa_dict_bykind[ikind] = bonds_aa_dict
                self._mass_aa_dict_bykind[ikind] = mass_aa_list

        # Create the aa_topology for each kind
        topology_aa_dict = defaultdict()
        for ikind in range(self._nkinds_aa):
            topology_aa_dict[ikind] = Topology()
            for iatom in range(self._natoms_aa[ikind]):
                topology_aa_dict[ikind].add_vertex(iatom)
            for ibond, values in self._bonds_aa_dict_bykind[ikind].items():
                for jbond in values:
                    topology_aa_dict[ikind].add_edge([ibond, jbond])

        # Get the distribution of molecule types in the top file.
        idx = indices_lines["molecules"][0] + 1
        while True:
            try:
                iline = topo_contents[idx]
                itype, nmols = iline.split()
                self._distribution_molecules.append([itype, int(nmols)])
                idx += 1
            except IndexError:
                break

        # Map from name of mols in aa model to ikind.
        ikind = 0
        for idx in indices_lines["moleculetype"]:
            while True:
                iline = topo_contents[idx + 1]
                if re.match(r"^;", iline):
                    idx += 1
                    continue
                try:
                    iname, tmp = iline.split()
                    idx += 1
                    self._map_aanamemol_to_ikind[iname] = ikind
                except ValueError:
                    break
            ikind += 1

        # Debug information
        if self._debug:
            for ikind in range(self._nkinds_aa):
                topology_aa_dict[ikind].draw_graph_forest_pygraphviz("topology_model{}aa".format(ikind))

        return topology_aa_dict

    # #########################################################################
    def _create_cg_topology(self):

        # Deduce CG bonds --> Find for the intersection of bonds within the aa representation
        # Get all aa bonds inside a bead
        aa_bonds_in_cg_bead_bykind = defaultdict()
        for ikind in range(self._nkinds_aa):
            aa_bonds_in_cg_bead = defaultdict(list)
            for ibead, iatoms_in_bead in enumerate(self._aggregate_index[ikind]):
                for iatom in iatoms_in_bead:
                    for jatom in self._bonds_aa_dict_bykind[ikind][iatom]:
                        b = [min(iatom, jatom), max(iatom, jatom)]
                        aa_bonds_in_cg_bead[ibead].append(b)
                # Remove duplicates
                aa_bonds_in_cg_bead[ibead] = list(set(tuple(sub) for sub in aa_bonds_in_cg_bead[ibead]))
            aa_bonds_in_cg_bead_bykind[ikind] = aa_bonds_in_cg_bead

        # Search for aa bonds between cg beads (intersection)
        for ikind in range(self._nkinds_aa):
            cg_bonds = []
            for ibead in range(self._nbeads_cg[ikind]):
                a = set(aa_bonds_in_cg_bead_bykind[ikind][ibead])
                for jbead in range(ibead+1, self._nbeads_cg[ikind]):
                    b = set(aa_bonds_in_cg_bead_bykind[ikind][jbead])
                    if len(a & b) > 0:
                        cg_bonds.append([ibead, jbead])
            self._cg_bond_bykind[ikind] = cg_bonds

        # Create CG topology
        topology_cg_dict = defaultdict()
        itable_bond = 1
        itable_ang = 1
        itable_dih = 1
        for ikind in range(self._nkinds_aa):
            topology_cg_dict[ikind] = Topology()
            for ibead in range(self._nbeads_cg[ikind]):
                topology_cg_dict[ikind].add_vertex(ibead)

            topology_cg_dict[ikind].set_name(self._type_bead_bykind[ikind])
            topology_cg_dict[ikind].set_type(self._type_name_bykind[ikind])

            bond_cg_type_dict = defaultdict(list)
            for item in self._cg_bond_bykind[ikind]:
                icg = item[0]
                jcg = item[1]
                topology_cg_dict[ikind].add_edge([icg, jcg])
                icgtype = self._type_bead_bykind[ikind][icg]
                jcgtype = self._type_bead_bykind[ikind][jcg]
                icgname = self._type_name_bykind[ikind][icg]
                jcgname = self._type_name_bykind[ikind][jcg]
                # A-B
                label_bond = "-".join(sorted([icgtype, jcgtype]))
                # A1 B2
                label_name = " ".join(sorted([icgname, jcgname]))
                bond_cg_type_dict[label_bond].append(label_name)
                if label_bond not in self._bond_cg_numtables_dict.keys():
                    self._bond_cg_numtables_dict[label_bond] = itable_bond
                    itable_bond += 1
            self._bond_cg_type_dict_bykind[ikind] = bond_cg_type_dict

            self._cg_angle_bykind[ikind] = topology_cg_dict[ikind].get_allbends()
            angle_cg_type_dict = defaultdict(list)
            for item in self._cg_angle_bykind[ikind]:
                icg = item[0]
                jcg = item[1]
                kcg = item[2]
                icgtype = self._type_bead_bykind[ikind][icg]
                jcgtype = self._type_bead_bykind[ikind][jcg]
                kcgtype = self._type_bead_bykind[ikind][kcg]
                icgname = self._type_name_bykind[ikind][icg]
                jcgname = self._type_name_bykind[ikind][jcg]
                kcgname = self._type_name_bykind[ikind][kcg]
                # A-B-C
                ll = sorted([icgtype, kcgtype])
                kk = sorted([icgname, kcgname])
                label_angle = ll[0]+"-"+jcgtype+"-"+ll[1]
                # A1 B2 B3
                label_name = kk[0]+" "+jcgname+" "+kk[1]
                angle_cg_type_dict[label_angle].append(label_name)
                if label_angle not in self._angle_cg_numtables_dict.keys():
                    self._angle_cg_numtables_dict[label_angle] = itable_ang
                    itable_ang += 1
            self._angle_cg_type_dict_bykind[ikind] = angle_cg_type_dict

            self._cg_dihedral_bykind[ikind] = topology_cg_dict[ikind].get_alldihedrals()
            dihedral_cg_type_dict = defaultdict(list)
            for item in self._cg_dihedral_bykind[ikind]:
                icg = item[0]
                jcg = item[1]
                kcg = item[2]
                lcg = item[3]
                icgtype = self._type_bead_bykind[ikind][icg]
                jcgtype = self._type_bead_bykind[ikind][jcg]
                kcgtype = self._type_bead_bykind[ikind][kcg]
                lcgtype = self._type_bead_bykind[ikind][lcg]
                icgname = self._type_name_bykind[ikind][icg]
                jcgname = self._type_name_bykind[ikind][jcg]
                kcgname = self._type_name_bykind[ikind][kcg]
                lcgname = self._type_name_bykind[ikind][lcg]
                # A-B-C-D
                label_dihedral = icgtype + "-" + jcgtype + "-" + kcgtype + "-" + lcgtype
                label_name = icgname + " " + jcgname + " " + kcgname + " " + lcgname
                label_dihedral_inv = lcgtype + "-" + kcgtype + "-" + jcgtype + "-" + icgtype
                label_name_inv = lcgname + " " + kcgname + " " + jcgname + " " + icgname
                if label_dihedral in dihedral_cg_type_dict.keys():
                    dihedral_cg_type_dict[label_dihedral].append(label_name)
                elif label_dihedral_inv in dihedral_cg_type_dict.keys():
                    dihedral_cg_type_dict[label_dihedral_inv].append(label_name_inv)
                    label_dihedral = label_dihedral_inv
                else:
                    dihedral_cg_type_dict[label_dihedral].append(label_name)

                if label_dihedral not in self._dihedral_cg_numtables_dict.keys():
                    self._dihedral_cg_numtables_dict[label_dihedral] = itable_dih
                    itable_dih += 1
            self._dihedral_cg_type_dict_bykind[ikind] = dihedral_cg_type_dict

    # #########################################################################
    def _check_parameters(self, label, contents=None, start_idx=None, end_idx=None):

        if label == "general":

            if self._title_system is None or len(self._title_system.split()) == 0:
                m = "\n\t\tERROR: title_system\n"
                m += "\t\t\t <GENERAL> ... </GENERAL> error\n"
                m += "First line in GENERAL: title_system cannot be empty"
                print(m) if self._logger is None else self._logger.error(m)
                exit()

            try:
                self._nkinds_aa = int(self._nkinds_aa)
            except ValueError:
                m = "\n\t\tERROR: nkinds_aa\n"
                m += "\t\t\t <GENERAL> ... </GENERAL> error\n"
                m += "Second line in GENERAL: nkinds_aa must be integer"
                print(m) if self._logger is None else self._logger.error(m)
                exit()

            if not os.path.isfile(self._gro_filename_aa):
                m = "\n\t\tERROR: gro_filename_aa\n"
                m += "\t\t\t <GENERAL> ... </GENERAL> error\n"
                m += "Third line in GENERAL: gro_filename_aa file must exist."
                print(m) if self._logger is None else self._logger.error(m)
                exit()

            if not os.path.isfile(self._top_filename_aa):
                m = "\n\t\tERROR: top_filename_aa\n"
                m += "\t\t\t <GENERAL> ... </GENERAL> error\n"
                m += "Fourth line in GENERAL: top_filename_aa file must exist."
                print(m) if self._logger is None else self._logger.error(m)
                exit()

        elif label == "names":

            if len(self._name_topol) != self._nkinds_aa or \
               len(self._name_cg) != self._nkinds_aa or \
               len(self._xml_filename_map_aa_cg) != self._nkinds_aa:
                m = "\n\t\tERROR: top_filename_aa\n"
                m += "\t\t\t <NAMES> ... </NAMES> error\n"
                m += "Section NAMES: The number of lines must be the same than the number of molecule kinds."
                print(m) if self._logger is None else self._logger.error(m)
                exit()

        elif label == "maps":
            if len(self._nbeads_cg) != self._nkinds_aa or \
               len(self._natoms_aa) != self._nkinds_aa or \
               len(self._nchains_aa_for_kind) != self._nkinds_aa:
                m = "\n\t\tERROR: top_filename_aa\n"
                m += "\t\t\t <MAPS> ... </MAPS> error\n"
                m += "Section MAPS: The number of lines must be the same than the number of molecule kinds."
                print(m) if self._logger is None else self._logger.error(m)
                exit()

            m = "\n\t\t Number of molecule/polymer types: {}\n".format(self._nkinds_aa)
            for ikind in range(self._nkinds_aa):
                try:
                    self._nbeads_cg[ikind] = int(self._nbeads_cg[ikind])
                    self._natoms_aa[ikind] = int(self._natoms_aa[ikind])
                    self._nchains_aa_for_kind[ikind] = int(self._nchains_aa_for_kind[ikind])
                except ValueError:
                    m = "\n\t\tERROR: top_filename_aa\n"
                    m += "\t\t\t <MAPS> ... </MAPS> error\n"
                    m += "Section MAPS format: natoms_aa(int) nbeads_cg(int) nchains(int) for each kind"
                    print(m) if self._logger is None else self._logger.error(m)
                    exit()

                m += "\t\t\t Type {0:3d}:    name={1:<8s} AA_atoms={2:<8d} CG_beads={3:<8d} Nmols={4:<6d}\n".\
                    format(ikind, self._name_cg[ikind], self._natoms_aa[ikind],
                           self._nbeads_cg[ikind], self._nchains_aa_for_kind[ikind])
            print(m) if self._logger is None else self._logger.info(m)

        elif label == "beads_previous":

            idx = start_idx + 1
            ikind_list = []
            while idx < end_idx:
                try:
                    ikind_list.append(int(contents[idx].split()[0]))
                    nlines = int(contents[idx].split()[1])
                    idx += nlines + 1
                except ValueError:
                    m = "\n\t\tERROR: \n"
                    m += "\t\t\t <BEADS> ... </BEADS> error\n"
                    m += "Section BEADS format header: idx_kind(int) number_cg_beads(int) for this kind\n"
                    m += "Last read: {}".format(contents[idx])
                    print(m) if self._logger is None else self._logger.error(m)
                    exit()

            if len(ikind_list) != len(set(ikind_list)):
                m = "\n\t\tERROR: \n"
                m += "\t\t\t <BEADS> ... </BEADS> error\n"
                m += "Section BEADS format: kinds must be different\n"
                m += "Last read: {}".format(contents[idx])
                print(m) if self._logger is None else self._logger.error(m)
                exit()

        elif label == "beads":
            pass

    # #########################################################################
    def write_map_cg_xml(self):
        """
        Write the map for each CG model to run CG simulations

        :return:
        """
        for ikind in range(self._nkinds_aa):
            with open(self._xml_filename_cg[ikind], 'w') as fout:
                fout.writelines("<cg_molecule>\n")

                line = "  <name>{0:s}</name>\n".format(self._name_cg[ikind])
                fout.writelines(line)

                # NAME OF THE MOLECULE IN THE TOP FILE (AA REPRESENTATION)
                line = "  <ident>{0:s}</ident>\n".format(self._name_cg[ikind])
                fout.writelines(line)

                # TOPOLOGY --> DEFINING THE CG BEADS
                fout.writelines("  <topology>\n")
                fout.writelines("   <cg_beads>\n")

                # WRITE CG BEADS
                # idx_global = 0
                for ibead in range(0, self._nbeads_cg[ikind]):
                    fout.writelines("     <cg_bead>\n")
                    fout.writelines("       <name>{}</name>\n".format(self._type_name_bykind[ikind][ibead]))
                    fout.writelines("       <type>{}</type>\n".format(self._type_bead_bykind[ikind][ibead]))
                    fout.writelines("       <mapping>{}</mapping>\n".format("UNITY"))

                    resid = "1 "
                    resname = self._name_cg[ikind]
                    atomname = self._type_name_bykind[ikind][ibead]
                    string_beads = resid.strip() + ":" + resname.strip() + ":" + atomname.strip()

                    fout.writelines("       <beads>{}</beads>\n".format(string_beads))
                    fout.writelines("     </cg_bead>\n")

                fout.writelines("   </cg_beads>\n")

                # BONDING OF THE CG MODEL
                fout.writelines("   <cg_bonded>\n")
                # BONDS
                for type_bond, bond_names in self._bond_cg_type_dict_bykind[ikind].items():
                    fout.writelines("     <bond>\n")
                    fout.writelines("       <name>{}</name>\n".format(type_bond))
                    fout.writelines("     <beads>\n")
                    for ibond in bond_names:
                        fout.writelines("       {}\n".format(ibond))
                    fout.writelines("     </beads>\n")
                    fout.writelines("     </bond>\n")

                # ANGLES
                for type_angle, angle_names in self._angle_cg_type_dict_bykind[ikind].items():
                    fout.writelines("     <angle>\n")
                    fout.writelines("       <name>{}</name>\n".format(type_angle))
                    fout.writelines("     <beads>\n")
                    for iangle in angle_names:
                        fout.writelines("       {}\n".format(iangle))
                    fout.writelines("     </beads>\n")
                    fout.writelines("     </angle>\n")

                # DIHEDRALS
                for type_dihedral, dihedral_names in self._dihedral_cg_type_dict_bykind[ikind].items():
                    fout.writelines("     <dihedral>\n")
                    fout.writelines("       <name>{}</name>\n".format(type_dihedral))
                    fout.writelines("     <beads>\n")
                    for idih in dihedral_names:
                        fout.writelines("       {}\n".format(idih))
                    fout.writelines("     </beads>\n")
                    fout.writelines("     </dihedral>\n")

                fout.writelines("   </cg_bonded>\n")
                fout.writelines("  </topology>\n")

                # MAPS
                fout.writelines("  <maps>\n")
                fout.writelines("    <map>\n")
                fout.writelines("      <name>{}</name>\n".format("UNITY"))
                fout.writelines("      <weights>1</weights>")
                fout.writelines("    </map>\n")
                fout.writelines("  </maps>\n")

                fout.writelines("</cg_molecule>\n")

    # #########################################################################
    def write_settings_initialdist(self, bond_min=0.1, bond_max=0.6, bond_step=0.02,
                                   angle_min=0.0, angle_max=3.1, angle_step=0.0348,
                                   dih_min=-3.1, dih_max=3.1, dih_step=0.08,
                                   nb_min=0.0, nb_max=2.0, nb_step=0.01, temp=500,
                                   path_gmxrc="/opt/votca_2022/bin/GMXRC.bash"):

        with open("settings_map.xml", 'w') as fout:

            lines = "<cg>\n\n"
            idx_bond = 1
            for ikind in range(self._nkinds_aa):
                for ibond in self._bond_cg_type_dict_bykind[ikind].keys():
                    lines += "  <bonded>\n"
                    lines += "      <name>{}</name>\n".format(ibond)
                    lines += "      <min>{}</min>\n".format(bond_min)
                    lines += "      <max>{}</max>\n".format(bond_max)
                    lines += "      <step>{}</step>\n".format(bond_step)
                    lines += "      <inverse>\n"
                    lines += "          <!-- Target distribution -->\n"
                    lines += "          <target>{}</target>\n".format("./"+ibond+".dist.tgt")
                    lines += "      </inverse>\n"
                    lines += "  </bonded>\n\n"
                    idx_bond += 1

            for ikind in range(self._nkinds_aa):
                for iangle in self._angle_cg_type_dict_bykind[ikind].keys():
                    lines += "  <bonded>\n"
                    lines += "      <name>{}</name>\n".format(iangle)
                    lines += "      <min>{}</min>\n".format(angle_min)
                    lines += "      <max>{}</max>\n".format(angle_max)
                    lines += "      <step>{}</step>\n".format(angle_step)
                    lines += "      <inverse>\n"
                    lines += "          <!-- Target distribution -->\n"
                    lines += "          <target>{}</target>\n".format("./"+iangle+".dist.tgt")
                    lines += "      </inverse>\n"
                    lines += "  </bonded>\n\n"

            for ikind in range(self._nkinds_aa):
                for idih in self._dihedral_cg_type_dict_bykind[ikind].keys():
                    lines += "  <bonded>\n"
                    lines += "      <name>{}</name>\n".format(idih)
                    lines += "      <min>{}</min>\n".format(dih_min)
                    lines += "      <max>{}</max>\n".format(dih_max)
                    lines += "      <step>{}</step>\n".format(dih_step)
                    lines += "      <inverse>\n"
                    lines += "          <!-- Target distribution -->\n"
                    lines += "          <target>{}</target>\n".format("./"+idih+".dist.tgt")
                    lines += "      </inverse>\n"
                    lines += "  </bonded>\n\n"

            label_cg_for_non_bonded = []
            for ikind in range(0, self._nkinds_aa):
                ll = self._cg_dict[ikind].keys()
                label_cg_for_non_bonded += ll

            for idx_ilabel in range(0, len(label_cg_for_non_bonded)):
                for idx_jlabel in range(idx_ilabel, len(label_cg_for_non_bonded)):
                    ilabel = label_cg_for_non_bonded[idx_ilabel]
                    jlabel = label_cg_for_non_bonded[idx_jlabel]
                    lines += "  <non-bonded>\n"
                    lines += "      <!-- name of the interaction -->\n"
                    lines += "      <name>{}{}</name>\n".format(ilabel, jlabel)
                    lines += "      <!-- types involved in this interaction -->\n"
                    lines += "      <type1>{}</type1>\n".format(ilabel)
                    lines += "      <type2>{}</type2>\n".format(jlabel)
                    lines += "      <!-- dimension + grid spacing of tables for calculations -->\n"
                    lines += "      <min>{}</min>\n".format(nb_min)
                    lines += "      <max>{}</max>\n".format(nb_max)
                    lines += "      <step>{}</step>\n".format(nb_step)
                    lines += "      <inverse>\n"
                    lines += "          <target>./{}{}.dist.tgt </target>\n".format(ilabel, jlabel)
                    lines += "      </inverse>\n"
                    lines += "  </non-bonded>\n\n"

            lines += "  <inverse>\n "
            lines += "      <scriptpath>./</scriptpath>\n"
            lines += "      <!-- {} * 0.00831451 gromacs units -->\n".format(temp)
            lines += "      <kBT>{}</kBT>\n".format(temp*0.00831451)
            lines += "      <dist_min>0.001</dist_min>\n"
            lines += "      <program>gromacs</program>\n"
            lines += "      <gromacs>\n"
            lines += "          <gmxrc> #PATH_TO_GMXRC# </gmxrc>\n"
            lines += "      </gromacs>\n"
            lines += "  </inverse>\n\n"

            lines += "</cg>\n\n"

            fout.writelines(lines)

    # #########################################################################
    def write_settings_cg(self, bond_min=0.1, bond_max=0.6, bond_step=0.02,
                          angle_min=0.0, angle_max=3.1, angle_step=0.0348,
                          dih_min=-3.1, dih_max=3.1, dih_step=0.08,
                          nb_min=0.0, nb_max=2.0, nb_step=0.01, temp=500,
                          path_gmxrc="/opt/votca_2022/bin/GMXRC.bash", max_iter=500):

        with open("settings_cg.xml", 'w') as fout:

            lines = "<cg>\n\n"
            idx_bond = 1
            for ikind in range(self._nkinds_aa):
                for ibond in self._bond_cg_type_dict_bykind[ikind].keys():
                    lines += "  <bonded>\n"
                    lines += "      <name>{}</name>\n".format(ibond)
                    lines += "      <min>{}</min>\n".format(bond_min)
                    lines += "      <max>{}</max>\n".format(bond_max)
                    lines += "      <step>{}</step>\n".format(bond_step)
                    lines += "      <inverse>\n"
                    lines += "          <!-- Target distribution -->\n"
                    lines += "          <target>{}</target>\n".format("./"+ibond+".dist.tgt")
                    lines += "          <gromacs>\n"
                    lines += "             <table>{}</table>\n".format("table_b"+str(idx_bond)+".xvg")
                    lines += "          </gromacs>\n"
                    lines += "      </inverse>\n"
                    lines += "  </bonded>\n\n"
                    idx_bond += 1

            idx_angle = 1
            for ikind in range(self._nkinds_aa):
                for iangle in self._angle_cg_type_dict_bykind[ikind].keys():
                    lines += "  <bonded>\n"
                    lines += "      <name>{}</name>\n".format(iangle)
                    lines += "      <min>{}</min>\n".format(angle_min)
                    lines += "      <max>{}</max>\n".format(angle_max)
                    lines += "      <step>{}</step>\n".format(angle_step)
                    lines += "      <inverse>\n"
                    lines += "          <!-- Target distribution -->\n"
                    lines += "          <target>{}</target>\n".format("./"+iangle+".dist.tgt")
                    lines += "          <gromacs>\n"
                    lines += "             <table>{}</table>\n".format("table_a"+str(idx_angle)+".xvg")
                    lines += "          </gromacs>\n"
                    lines += "      </inverse>\n"
                    lines += "  </bonded>\n\n"
                    idx_angle += 1

            idx_dihedral = 1
            for ikind in range(self._nkinds_aa):
                for idih in self._dihedral_cg_type_dict_bykind[ikind].keys():
                    lines += "  <bonded>\n"
                    lines += "      <name>{}</name>\n".format(idih)
                    lines += "      <min>{}</min>\n".format(dih_min)
                    lines += "      <max>{}</max>\n".format(dih_max)
                    lines += "      <step>{}</step>\n".format(dih_step)
                    lines += "      <inverse>\n"
                    lines += "          <!-- Target distribution -->\n"
                    lines += "          <target>{}</target>\n".format("./"+idih+".dist.tgt")
                    lines += "          <gromacs>\n"
                    lines += "             <table>{}</table>\n".format("table_d"+str(idx_dihedral)+".xvg")
                    lines += "          </gromacs>\n"
                    lines += "      </inverse>\n"
                    lines += "  </bonded>\n\n"
                    idx_dihedral += 1

            label_cg_for_non_bonded = []
            do_potential_list = []
            for ikind in range(0, self._nkinds_aa):
                ll = self._cg_dict[ikind].keys()
                label_cg_for_non_bonded += ll

            for idx_ilabel in range(0, len(label_cg_for_non_bonded)):
                for idx_jlabel in range(idx_ilabel, len(label_cg_for_non_bonded)):
                    do_potential_list.append(0)

            idx_nb = 0
            for idx_ilabel in range(0, len(label_cg_for_non_bonded)):
                for idx_jlabel in range(idx_ilabel, len(label_cg_for_non_bonded)):
                    ilabel = label_cg_for_non_bonded[idx_ilabel]
                    jlabel = label_cg_for_non_bonded[idx_jlabel]
                    lines += "  <non-bonded>\n"
                    lines += "      <!-- name of the interaction -->\n"
                    lines += "      <name>{}{}</name>\n".format(ilabel, jlabel)
                    lines += "      <!-- types involved in this interaction -->\n"
                    lines += "      <type1>{}</type1>\n".format(ilabel)
                    lines += "      <type2>{}</type2>\n".format(jlabel)
                    lines += "      <!-- dimension + grid spacing of tables for calculations -->\n"
                    lines += "      <min>{}</min>\n".format(nb_min)
                    lines += "      <max>{}</max>\n".format(nb_max)
                    lines += "      <step>{}</step>\n".format(nb_step)
                    lines += "      <inverse>\n"
                    lines += "          <target>./{}{}.dist.tgt </target>\n".format(ilabel, jlabel)
                    lines += "          <do_potential>"
                    for i in range(len(do_potential_list)):
                        if i == idx_nb:
                            lines += " {}".format(str(1))
                        else:
                            lines += " {}".format(str(do_potential_list[i]))
                    lines += "</do_potential>\n"
                    lines += "          <post_update></post_update>\n"
                    lines += "          <!-- additional post processing of U after dU added to potential -->\n"
                    lines += "          <post_add>convergence</post_add>\n"
                    lines += "          <gromacs>\n"
                    lines += "              <table>table_{}_{}.xvg</table>\n".format(ilabel, jlabel)
                    lines += "          </gromacs>\n"
                    lines += "      </inverse>\n"
                    lines += "  </non-bonded>\n\n"
                    idx_nb += 1

            lines += "  <inverse>\n "
            lines += "      <scriptpath>./</scriptpath>\n"
            lines += "      <!-- {} * 0.00831451 gromacs units -->\n".format(temp)
            lines += "      <kBT>{}</kBT>\n".format(temp*0.00831451)
            lines += "      <dist_min>0.001</dist_min>\n"
            lines += "      <program>gromacs</program>\n"
            lines += "      <gromacs>\n"
            lines += "          <conf> #GRO_CG_FILENAME# </conf>\n"
            lines += "          <topol_in> #TOPO_CG_FILENAME# </topol_in>\n"
            lines += "          <index> #INDEX_CG_FILENAME# </index>\n"
            lines += "          <gmxrc> #PATH_TO_GMXRC# </gmxrc>\n"
            lines += "          <equi_time>20</equi_time>\n"
            lines += "          <table_bins>0.002</table_bins>\n"
            lines += "          <pot_max>1000000</pot_max>\n"
            lines += "          <table_end>3.0</table_end>\n"
            lines += "          <grompp>\n"
            lines += "             <bin> #PATH_TO_GMX_GROMPP# </bin>\n"
            lines += "             <opts>-maxwarn 6</opts>\n"
            lines += "          </grompp>\n"
            lines += "          <mdrun>\n"
            lines += "             <command> #PATH_TO_GMX_MDRUN# </command>\n"
            lines += "             <opts>-ntmpi 48 -rdd 1.0</opts>\n"
            lines += "          </mdrun>\n"
            lines += "          <g_energy>\n"
            lines += "             <bin> #PATH_TO_GMX_ENERGY# </bin>\n"
            lines += "          </g_energy>\n"
            lines += "      </gromacs>\n"
            lines += "      <map> #XML_CG_FILE# </map>\n"
            lines += "      <filelist>  #GRO_CG_FILENAME# grompp.mdp #TOPO_CG_FILENAME# table.xvg " \
                     "#INDEX_CG_FILENAME#  </filelist>\n"
            lines += "      <iterations_max>{}</iterations_max>".format(max_iter)
            lines += "      <convergence_check>\n"
            lines += "          <type>default</type>\n"
            lines += "          <limit>0.001</limit>\n"
            lines += "      </convergence_check>\n"
            lines += "      <!-- ibm: inverse boltzmann imc: inverse monte carlo -->\n"
            lines += "      <method>ibi</method>\n"
            lines += "      <!-- write log to this file -->\n"
            lines += "      <log_file>inverse.log</log_file>\n"
            lines += "      <!-- write restart step to this file -->\n"
            lines += "      <restart_file>restart_points.log</restart_file>\n"
            lines += "  </inverse>\n\n"
            lines += "</cg>\n\n"

            fout.writelines(lines)
