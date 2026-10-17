# Run from the amp-atlas project directory in an existing PyMOL installation.
# Experimental PDB 2MAG, first model, chain A. Cite the original PDB study.
load data/structures/2MAG_41.pdb, magainin
hide everything
show cartoon, magainin
show sticks, magainin and not backbone
spectrum count, rainbow, magainin
bg_color white
set cartoon_fancy_helices, 1
set stick_radius, 0.15
set ray_opaque_background, off
orient magainin
# Uncomment to render a publication-style transparent PNG:
# ray 1800, 1200
# png magainin_2MAG.png, dpi=300
