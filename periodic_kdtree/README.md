# Fibronectins

This is a coarse-grained model of a fibronectin monomer. The monomer undergoes unfolding thus exposing hydrophobic binding sites that can interact with other protein monomers and thus drive the self-assembly process. The resulting structure is a fibrillar matrix, its morphology depending on various properties that will be investigated in this project.

---

In order to run this project, download and make a local LAMMPS executable.

---

Make the system as:

```python 
python fibronectin.py
```

This will create the topology file and the input files that can be fed directly to the LAMMSP executable.
 
Run it with:

```bash
./lmp_serial < in.fibro
```

---
In order to visualise the trajectory output, install the VMD package and do:

```bash
vmd fibro.xyz -e vmdfibro
```

