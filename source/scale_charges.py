from metaltool import Forcefield, get_atomnames

def main(forcefield_file:str, files_for_scaling:list, scale_factor:float):
    """Generates a new forcefield file scaling the charges of some atoms. The new forcefield will not have
       atomic polarizabilities.

    Args:
        forcefield_file (str): Path to the file containing the forcefield
        files_for_scaling (list): List of paths with the structures whose atoms' charges are to be scaled
        scale_factor (float): Multiplicative scale factor
    """

    forcefield = Forcefield.from_file(forcefield_file)
    names_for_scaling = []

    for file in files_for_scaling:

        names_for_scaling += get_atomnames(file)
    
    ## Remove possible duplicates
    names_for_scaling = set(names_for_scaling)

    ## Check if all atoms are in the forcefield
    names_in_ff = set([a.name for a in forcefield.atoms])
    for name_s in names_for_scaling:
        assert name_s in names_in_ff, f"Atom {name_s} is not in the forcefield"

    ## Perform the scaling

    for atom in forcefield.atoms:
        if atom.name in names_for_scaling:
            current_charge = float(atom.charge)
            new_charge = current_charge * scale_factor
            atom.charge = str(new_charge)

    ## Remove the polarization information from the forcefield
            
    for atom in forcefield.atoms:
        atom.polarizability = ""

    ## Save the new forcefield
            
    forcefield.write_file(f"scaled_{forcefield_file}",
                          comment=f"## Written with metaltool\n## The charge of atoms {' '.join(names_for_scaling)} was scaled\n## The scaling factor was {scale_factor}")

import sys

if __name__ == '__main__':

    scale = float(sys.argv[0])
    ff_file = sys.argv[1]
    structure_files = sys.argv[2:]

    main(ff_file, structure_files, scale)







    