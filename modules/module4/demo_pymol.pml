# BIOT 6900 demo: acetazolamide (3HS4) in the AlphaFold model of carbonic anhydrase II
load structures/AF-P00918_aligned_to_3HS4.pdb, af
load structures/3HS4.pdb, crystal
hide everything
show cartoon, af
color orange, af and b < 50
color yellow, af and b > 50 and b < 70
color lightblue, af and b > 70 and b < 90
color blue, af and b > 90
show sticks, crystal and resn AZM and resi 701
color magenta, crystal and resn AZM and elem C
show spheres, crystal and resn ZN
set sphere_scale, 0.4
show surface, af within 8 of (crystal and resn AZM and resi 701)
set transparency, 0.45
orient crystal and resn AZM and resi 701
zoom crystal and resn AZM and resi 701, 8
