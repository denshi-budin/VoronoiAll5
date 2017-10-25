# 3DVoronoiTess-5.0
Creating 3D Voronoi tessellations specimens for molecular dynamic simulation.

Work with Windows(cygwin), Mac and Linux (ubuntu).

## Getting Started

* To compile just type:
```
make voronoi
```
* Edit Configuration File (Configuration.cfg):
```
Lattice : (Diamond, Fcc, Hcp, Bcc)
Scale: (lattice scale in Angstrom)
Element: (Element symbol. e.g: Si, Cu, etc)
Mass: (Element mass)
```

* Running :
```
./voronoi [(optianal) Random Seed]
```

## Output files 

Output files can be found in output folder:

~in.mddata ; specimen's data format for LAMMPS

sample.xyz ; specimen's data format in xyz

## Authors

* **Ahmad Ehsan Mohd Tamidi** - *Initial work*
