"""
This script changes atom names in the gro and top file.
In order to do the map from AA to CG, the atoms in the AA
are identified as:
   nres:nameres:nameat

In this example the map is:
  3 atoms AA --> 1 bead centered in the com of these atoms.
  According to the input files

The first bead should be:
       1:PES:_CH3 1:PES:_CH2 2:PES:_CH2
The second and following are:
       2:PES:_CH2 3:PES:_CH2 3:PES:_CH2  Here there is a problem
       4:PES:_CH2 4:PES:_CH2 5:PES:_CH2  Here there is a problem
       5:PES:_CH2 6:PES:_CH2 6:PES:_CH2  Here there is a problem
and the last is defined as:
       33:PES:_CH2 34:PES:_CH2 34:PES:_CH3

To avoid problems we give a different name to each atom in the residue.
"""
import re
from collections import defaultdict

# ============= INPUTS ========================
faa_name_gro = "confout.part0001.gro"
faa_name_top = "PE-iPP_51mon_20-20Ch_residues_noH_replicate.top"
faa_name_gro_new = "conf_label.gro"
faa_name_top_new = "conf_label.top"
# ============= INPUTS ========================

letter = ["A", "B", "C", "D", "E", "F", "G", "H"]

# Format topopology ATOM
fmt_top = {'ATOM': "{idx:>6d}    {attype:7s}{resid:>7d}"
                   "{restypename:>7s}{atomname:>7s}"
                   "{chargegroup:>7d}{charge:>11.8f} "
                   "{mass:>11.8f}   {comment:s}\n"}

# Read topology file
with open(faa_name_top, 'r') as faa_top:

    contents_top = faa_top.readlines()

    ipos_s_atoms_labels = list()
    ipos_e_atoms_labels = list()
    residues_list = list()
    ipos_s_moleculetype_labels = list()
    ipos_s_molecules_labels = list()

    # Find position of the [atoms], [ moleculetype ] and [ molecules ] labels in the topology file ==============
    for idx in range(len(contents_top)):
        if contents_top[idx] == "[ atoms ]\n":
            ipos_s_atoms_labels.append(idx)
        if contents_top[idx] == "[ moleculetype ]\n":
            ipos_s_moleculetype_labels.append(idx)
        if contents_top[idx] == "[ molecules ]\n":
            ipos_s_molecules_labels.append(idx)

    # Check that the number of [ moleculetype ] labels are the same that the nummber of [ atoms ] labels
    nmol_types = len(ipos_s_atoms_labels)
    if nmol_types != len(ipos_s_moleculetype_labels):
        print("ERROR: The number of [ atoms ] and [ moleculetype ] labels must be equal!!!!!.")
        print("        # atoms        : {}".format(len(ipos_s_atoms_labels)))
        print("        # moleculetype : {}".format(len(ipos_s_moleculetype_labels)))
        exit()

    # Take lines for atoms =============================
    atoms_lines = defaultdict(list)
    for itype_mol in range(0, nmol_types):
        idx_line = ipos_s_atoms_labels[itype_mol] + 1
        while True:

            iline = str(contents_top[idx_line])

            if re.match(r'^(\s*)\n', iline):
                ipos_e_atoms_labels.append(idx_line)
                break

            else:
                atoms_lines[itype_mol].append(iline)
                idx_line += 1

    # Name of each type
    names_by_mol = list()
    map_topname_to_labelresname = defaultdict()
    itype_mol = 0
    for idx_moltype in ipos_s_moleculetype_labels:
        idx_tmp = idx_moltype
        while True:
            try:
                iline = contents_top[idx_tmp].split()[0]
            except IndexError:
                names_by_mol.append(iline)
                map_topname_to_labelresname[iline] = letter[itype_mol]
                itype_mol += 1
                break
            idx_tmp += 1

    # Natoms for each molecule
    natoms_by_type = defaultdict(list)
    for itype_mol in range(0, nmol_types):
        natoms_by_type[names_by_mol[itype_mol]] = 0


# Read the gro
# Open topology file
with open(faa_name_gro, 'r') as faa_gro:
    contents_gro = faa_gro.readlines()

# Write topology file =====================
with open(faa_name_top_new, 'w') as faa_top_new:

    # Write lines before first  [moleculetype] (topology)
    for i in contents_top[:ipos_s_moleculetype_labels[0]]:
        faa_top_new.writelines(i)

    resid_count_by_mol = defaultdict()
    for itype_mol in range(0, nmol_types):
        namemol = names_by_mol[itype_mol]
        # Write lines from [ moleculetype ] to [ atoms ]
        lrange = ipos_s_moleculetype_labels[itype_mol]
        rrange = ipos_s_atoms_labels[itype_mol]+1
        if itype_mol == 0:
            for i in contents_top[lrange:rrange]:
                faa_top_new.writelines(i)
        else:
            faa_top_new.writelines("[ atoms ]")

        # Modify and write the resname of each line in Atoms
        resid_count = 1
        atlocal_count = 0
        for iline in atoms_lines[itype_mol]:
            fields = iline.split()
            if fields[0].find(";") != -1:
                # It is a comment --> write line
                faa_top_new.writelines(iline)
                continue
            if resid_count != int(fields[2]):
                atlocal_count = 1
                resid_count = int(fields[2])
            else:
                atlocal_count += 1
            fields[4] = "{0:1s}{1:04d}".format(map_topname_to_labelresname[namemol], atlocal_count)

            tmp = ' '.join(fields[8:])

            jline = fmt_top['ATOM'].format(idx=int(fields[0]),
                                           attype=fields[1],
                                           resid=int(fields[2]),
                                           restypename=fields[3],
                                           atomname=fields[4],
                                           chargegroup=int(fields[5]),
                                           charge=float(fields[6]),
                                           mass=float(fields[7]),
                                           comment=tmp
                                            )
            # Write new lines
            faa_top_new.writelines(jline)
            natoms_by_type[names_by_mol[itype_mol]] += 1

        resid_count_by_mol[names_by_mol[itype_mol]] = resid_count
        faa_top_new.writelines("\n")

        if itype_mol != nmol_types - 1:
            lrange = ipos_e_atoms_labels[itype_mol]
            rrange = ipos_s_atoms_labels[itype_mol+1]
            for i in contents_top[lrange:rrange]:
                faa_top_new.writelines(i)
        else:
            lrange = ipos_e_atoms_labels[itype_mol]
            for i in contents_top[lrange:]:
                faa_top_new.writelines(i)

    # Indexes for each molecule
    idx_molecules = ipos_s_molecules_labels[0] + 1
    idx_atom = 0
    idxatom_to_labelCG = defaultdict()
    while True:
        try:
            iline = contents_top[idx_molecules]
            inamemol, inmols = iline.split()
            for imol in range(0, int(inmols)):
                for iat in range(natoms_by_type[inamemol]):
                    idxatom_to_labelCG[idx_atom] = map_topname_to_labelresname[inamemol]
                    idx_atom += 1
        except IndexError:
            break
        idx_molecules += 1


# Write coordinates files =====================
with open(faa_name_gro_new, 'w') as faa_gro_new:

    fmt_gro = {'ATOM': "{resid:>5d}{restypename:5s}{atomname:>5s}"
                       "{idx:>5d}{xcoord:>8.3f}{ycoord:>8.3f}{zcoord:8.3f}\n"}

    # Write two first lines (coordinates)
    for i in contents_gro[0:2]:
        faa_gro_new.writelines(i)

    # Modify and write the resname of each line in Atoms (coordinates)
    resid_count  = 1
    atlocal_count = 0
    for iline in contents_gro[2:-1]:

        if resid_count != int(iline[0:5]):
            atlocal_count = 1
            resid_count = int(iline[0:5])
        else:
            atlocal_count += 1
        atomname = "{0:1s}{1:04d}".format(idxatom_to_labelCG[int(iline[15:20])-1], atlocal_count)
        jline = fmt_gro['ATOM'].format(resid=int(iline[0:5]),
                                       restypename=iline[5:10],
                                       atomname=atomname,
                                       idx=int(iline[15:20]),
                                       xcoord=float(iline[20:28]),
                                       ycoord=float(iline[28:36]),
                                       zcoord=float(iline[36:44]))
        faa_gro_new.writelines(jline)

    # Write box line (coordinates)
    faa_gro_new.writelines(contents_gro[-1])

with open("info.log", 'w') as flog:

    print("Number of systems: {}".format(nmol_types))
    i = 0
    for ikey, ivalue in natoms_by_type.items():
        print("\t Mol {0:<4d} --> name: {1:<10s} natoms: {2:<8d} nresidues: {3:<8d}".format(i, ikey,
                                                                                         ivalue,
                                                                                         resid_count_by_mol[ikey]))
        i += 1

    flog.writelines("Number of systems: {}\n".format(nmol_types))
    i = 0
    for ikey, ivalue in natoms_by_type.items():
        flog.writelines("\t Mol {0:<4d} --> name: {1:<10s} natoms: {2:<8d} nresidues: {3:<8d}\n".format(i, ikey,
                                                                                                   ivalue,
                                                                                            resid_count_by_mol[ikey]))
        i += 1