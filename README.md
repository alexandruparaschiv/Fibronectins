[![Generic badge](https://img.shields.io/badge/status-in progress-green.svg)](https://shields.io/)
# Fibronectins


<p align="center">
  <img src="snapshots/Fn_mesh.png" width="532" height="400" title="hover text">
</p>

Fibronectins are structural proteins that provide a template for the growth of the extracellular matrix. The fibronectin protein is essentially a mechanochemical switch that changes conformation upon binding to cellular receptors. The monomer experiences cytoskeletal forces and  undergoes unfolding thus exposing hydrophobic binding sites. These can interact with other protein monomers and thus drive a self-assembly process. The resulting structure is a fibrillar matrix, its morphology depending on various properties that will be investigated in this project.

---

In order to run this project, download and make a local LAMMPS executable.

---

Make the system as:

```python 
python fibronectin.py
```

This will create the topology file and the input files that can be fed directly to the LAMMPS executable.
 
Run it with:

```bash
./lmp_serial < in.fibro
```

---
In order to visualise the trajectory output, install the VMD package and do:

```bash
vmd fibro.xyz -e vmdfibro
```
Alternatively, visualisation programs such as OVITO can be used.
