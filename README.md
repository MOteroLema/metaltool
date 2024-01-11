# metaltool
Python scripts to be used along fftool (https://github.com/paduagroup/fftool). The objective is to generate reduced topologies and forcefield parameters from the broad il.ff file provided in CL&P (https://github.com/paduagroup/clandp). The programs algo assigns a different name to every atom, to ensure compatibility with the METALWALLS software (https://gitlab.com/ampere2/metalwalls).

## Generating simplified topologies

Starting from a forcefield that contains many atoms, bonds, angles, dihedrals and improper dihedrals, and a set of structure files containing various molecules, it is possible to generate a reduced version of the forcefield that only contains information about the atoms present in the provided structure files. To do so, just run

    python metaltool.py forcefield.ff structure1.xyz structure2.zmat ... structureN.xyz

The supported formats for the structure files are xyz and zmat, and both can be used at the same time. Upon execution, the program will write a subset of the forcefield in a file called "new_forcefield.ff". Moreover, it will write new structure files called "new_structure.zmat" or "new_structure.xyz". In these new forcefield and structure files, no two atoms will have the same name, to ensure compatibility with the METALWALLS software, which requires different names for each atom in the given molecules. The resulting forcefield file will also have the value of the atomic polarization as the last item in each line in the ATOMS section.

## Scaling the charges in a forcefield

Charge scaling is a common practice in the simulation of ionic liquids, and serves as a way to mimick polarization effects in a computationally cheap way, without the need of Drude particles of induced dipoles. The ```scale_charges.py``` executable is a tool to generate a forcefield with scaled charges from an unscaled one. To generate the new forcefield, run the following code

    python scale_charges.py scale_factor forcefield.ff structure1.xyz ... structureN.zmat

This will generate a "scaled_forcefield.ff" file, where the charges of the atoms in the provided structures will be scaled by scale_factor (this number usually ranges between 0.7 and 0.8 for ionic liquids). Unlike ```metaltool.py```, this code won't eliminate any atoms from the forcefield nor rename them, in case only some molecules are desired to have their charges scaled. Moreover, the resulting forcefield will not have any information about atomic polarizability (it won't be written, not just set to 0.0), since charge scaling is meant to account for polarization effects.
