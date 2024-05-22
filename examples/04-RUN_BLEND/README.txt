# INPUT FILES
cd 01-PREPROCESS_GROFILES/
cp ../../../04-PE_AA/03-PE-iPP_blends/2022-PE-iPP_51Mon_160-160Ch_TRAPPETOX/06-NPT_CRes_200ns_500K/confout.part0001.gro .
cp ../../../04-PE_AA/03-PE-iPP_blends/2022-PE-iPP_51Mon_160-160Ch_TRAPPETOX/02-REPLICATE_TYPING_2_2_2/PE-iPP_51mon_20-20Ch_residues_noH_replicate.top .

# PREPROCESS ATOM NAMES WITH THE PYTHON SCRIPT AND CREATE A TPR
cd 01-PREPROCESS_GROFILES
python 00-script_AA_renameatoms_v02.py
gmx grompp -f minim.mdp -c conf_label.gro -p conf_label.top --maxwarn 2

# SETUP THE MAP IN THE PYTHON SCRIPT (Use topology library)
mkdir 02-MAP-iPP
cd ../02-MAP-iPP
source ~/Programacion/sandboxes/sandbox_common/bin/activate
python iPPCG3.py


