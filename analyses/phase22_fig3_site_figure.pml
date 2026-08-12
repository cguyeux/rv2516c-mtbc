# Figure 3 -- site de liaison a l'ADN candidat de l'unite B (P2.4/P2.6, modele Boltz-2 J1, sauvage).
# Chaine A = proteine (unite B, 60 residus, numerotation de fragment = H37Rv - 87).
# Chaines B/C = duplex d'ADN de 22 pb (operateur de BldC, PDB 6AMA).
# R113/H115/R118 (H37Rv) = resi 26/28/31 (fragment) -- offset verifie phase18_p2_6_interface_contacts.py.
# Rendu headless, depuis la racine du depot : pymol -cq analyses/phase22_fig3_site_figure.pml
# (chemins relatifs a la racine du depot ; PyMOL herite du repertoire de travail courant)

load résultats/p2_4_boltz/out_j1_unitB_dna/boltz_results_j1_unitB_dna/predictions/j1_unitB_dna/j1_unitB_dna_model_0.cif, j1

bg_color white
set ray_opaque_background, 0
set ray_shadows, 0
set antialias, 2
set cartoon_fancy_helices, 1
set valence, 0
hide everything

# --- ADN (chaines B/C) ---
show cartoon, j1 and chain B+C and polymer
set cartoon_ring_mode, 3
set cartoon_ring_finder, 1
set cartoon_nucleic_acid_mode, 4
color grey60, j1 and chain B+C

# --- proteine (unite B) ---
show cartoon, j1 and chain A and polymer
color skyblue, j1 and chain A
set cartoon_transparency, 0.25, j1 and chain A

# --- residus du site candidat ---
select site, j1 and chain A and resi 26+28+31
show sticks, site
set stick_radius, 0.22, site
color orange, site and elem C
util.cnc site

set label_size, 22
set label_color, black
set label_outline_color, white
set label_font_id, 7

# ======== Vue 1 : vue d'ensemble du complexe proteine-ADN (SANS labels, ambigus a cette echelle) ========
orient j1
ray 2200, 1800
png article/figures/figure3a_complex_overview.png, dpi=600

# ======== Vue 2 : gros plan sur le site de contact, labels decales pour ne pas se chevaucher ========
orient site
zoom site, 8.5
turn x, -10
set label_position, (-3, -4, 2), j1 and chain A and resi 26
label j1 and chain A and resi 26 and name CA, "Arg113"
set label_position, (7, 2, 2), j1 and chain A and resi 28
label j1 and chain A and resi 28 and name CA, "His115"
set label_position, (0, 6, 2), j1 and chain A and resi 31
label j1 and chain A and resi 31 and name CA, "Arg118"
ray 2200, 1800
png article/figures/figure3b_site_closeup.png, dpi=600
