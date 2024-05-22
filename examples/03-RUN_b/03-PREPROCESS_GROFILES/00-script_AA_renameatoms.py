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
faa_name_gro = "../02-REPLICATE_1_1_1/PE10-iPP10m_block_2Ch_residues_noH_replicate.gro"
faa_name_top = "../02-REPLICATE_1_1_1/PE10-iPP10m_block_2Ch_residues_noH_replicate.top"
faa_name_gro_new = "conf_label.gro"
faa_name_top_new = "conf_label.top"
nkinds = 1
# ============= INPUTS ========================

# Open topology file
with open(faa_name_top, 'r') as faa_top:

    contents_top = faa_top.readlines()

    ipos_atoms_labels = []
    ipos_start_labels = [0]
    atoms_lines = defaultdict(list)
    residues_list = []

    # Find position of the [atoms] label in the topology file ==============
    atoms_found = False
    for idx in range(len(contents_top)):
        if atoms_found:
            if len(contents_top[idx].strip()) < 2:
                ipos_start_labels.append(idx)
                atoms_found = False
        if contents_top[idx] == "[ atoms ]\n":
            ipos_atoms_labels.append(idx)
            atoms_found = True

    if len(ipos_atoms_labels) != nkinds:
        print("ERROR: Atoms must be appear only once in the topology file!!!!!.")
        print(ipos_atoms_labels)
        exit()

    # Take lines for atoms =============================
    for ikind in range(0, nkinds):

        idx_line = ipos_atoms_labels[ikind] + 1
        while True:

            iline = str(contents_top[idx_line])

            if re.match(r'^(\s*)\n', iline):
                ipos_atoms_end_label = idx_line
                break

            if re.match(r"^(\s*);", iline):
                idx_line += 1
            else:
                atoms_lines[ikind].append(iline)
                idx_line += 1

# Read the gro
# Open topology file
with open(faa_name_gro, 'r') as faa_gro:

    contents_gro = faa_gro.readlines()

# Write topology file =====================
with open(faa_name_top_new, 'w') as faa_top_new:

    fmt_top = {'ATOM': "{idx:>6d}    {attype:7s}{resid:>7d}"
                       "{restypename:>7s}{atomname:>7s}"
                       "{chargegroup:>7d}{charge:>11.8f} "
                       "{mass:>11.8f}   {comment:s}\n"}

    for ikind in range(nkinds):

        # Write lines before atoms (topology)
        for i in contents_top[ipos_start_labels[ikind]:ipos_atoms_labels[ikind] + 1]:
            faa_top_new.writelines(i)

        # Modify and write the resname of each line in Atoms
        resid_count = 1
        atlocal_count = 0
        for iline in atoms_lines[ikind]:
            fields = iline.split()

            if resid_count != int(fields[2]):
                atlocal_count = 1
                resid_count = int(fields[2])
            else:
                atlocal_count += 1
            fields[4] = "A{0:04d}".format(atlocal_count)

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


    for i in contents_top[ipos_atoms_end_label:]:
        faa_top_new.writelines(i)

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
        atomname = "A{0:04d}".format(atlocal_count)

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
