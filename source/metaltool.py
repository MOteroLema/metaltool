import sys
from collections import Counter
import os


class Atom():

    def __init__(self, name:str, type:str, mass:str, charge:str,
                 potential_type:str, potential_params:list, polarizability:str="0.0"):
        """Class to store atomic information

        Args:
            name (str): Name of the atom
            type (str): Type of the atom (this is what later defines which bond/angle/dihedral parameters are applied to it)
            mass (str): Mass of the atom (u)
            charge (str): Charge of the atom (e)
            potential_type (str): Either "lj" (Lennard-Jones) of "ft" (Fumi-Tosi)
            potential_params (list): List with the potential's parameters. For "lj" atoms, the list contains [sigma (A), epsilon (kJ/mol)], for 
                                     for "ft", the list contains [B (A^-1), A (kj/mol), C6 (kJ/(mol*A^-6)), C8 (kJ/mol*A^-8)]
            polarizability (str, optional): Polarizability of the atom (A^-3). Defaults to "0.0". If polarizabilities should not be written, input ""

        Raises:
            ValueError: If a non suported potential_type is provided
        """

        self.name = name
        self.atomtype = type
        self.mass = mass
        self.charge = charge
        self.polarizability = polarizability
        self._potential_params = potential_params

        if potential_type == "lj":
            assert len(potential_params)==2, f"Atom {name} of type {type} has been asigned a {potential_type} potential, which requires 2 parameters, and {len(potential_params)} were provided"
            sigma, epsilon = potential_params
            self.pottype = potential_type
            self.sigma = sigma
            self.epsilon = epsilon
        elif potential_type == "ft":
            assert len(potential_params)==4, f"Atom {name} of type {type} has been asigned a {potential_type} potential, which requires 4 parameters, and {len(potential_params)} were provided"
            self.B, self.A, self.C6, self.C8 = potential_params
            self.pottype = potential_type
        else:
            raise ValueError(f"Atom {name} of type {type} has been assigned a {potential_type} potential, which is not supported. Supported types are Lennard-Jones (lj) and Fumi-Tosi (ft)")
    
    @classmethod
    def parse_line(cls, line:str):
        """Defines an object from a string of characteristics

        Args:
            line (str): Line describing an atom, following the format of CL&P, with the option to have an extra parameter
                        parameter at the end of the line which will be the polarizability

        Raises:
            ValueError: If a potential type not supported by the Atom class is read from the line

        Returns:
            Atom: Object of the Atom class
        """

        items = line.strip().split()

        if items[4] == "lj":

            assert len(items)==7 or len(items)==8, f"Detected {len(items)} items in the line. Atoms of type lj only have either 7 or 8 items (without or with polarization)"
            name = items[0]
            atomtype = items[1]
            mass = items[2]
            charge = items[3]
            pottype = items[4]
            potential_params = [items[5], items[6]]

            ## Not all files will have polarization
            try:
                pol = items[7]
            except:
                pol = "0.0"

            return cls(name, atomtype, mass, charge, pottype, potential_params, pol)
        
        elif items[4] == "ft":

            assert len(items)==11 or len(items)==12, f"Detected {len(items)} items in the line. Atoms of type ft only have either 11 or 12 items (without or with polarization)"
            name = items[0]
            atomtype = items[1]
            mass = items[2]
            charge = items[3]
            pottype = items[4]
            potential_params = [items[5], items[6], items[7], items[8]]

            ## Not all files will have polarization
            try:
                pol = items[11]
            except:
                pol = 0

            return cls(name, atomtype, mass, charge, pottype, potential_params, pol)
        
        else:
            raise ValueError(f"This atom has been assigned a {items[4]} potential, which is not supported. Supported types are Lennard-Jones (lj) and Fumi-Tosi (ft)")

    def generate_line(self)->str:
        """Function to create lines for the forcefield

        Raises:
            ValueError: If the current potential type is not supported

        Returns:
            str: Line describing the atom, in the CL&P format
        """

        if self.pottype == "lj":

            contents = [self.name, self.atomtype, self.mass, self.charge, self.pottype, self.sigma, self.epsilon, self.polarizability]

        elif self.pottype == "ft":

            contents = [self.name, self.atomtype, self.mass, self.charge, self.pottype, self.B, self.A, self.C6, self.C8, self.B, self.B, self.polarizability]

        else:
            raise ValueError(f"Atom {self.name} of type {self.atomtype} has been assigned a {self.pottype} potential, which is not supported. Supported types are Lennard-Jones (lj) and Fumi-Tosi (ft)")

        line = "    ".join(contents)

        return line
    

class Bond():

    def __init__(self, type1:str, type2:str, bondtype:str, r_eq:str, k:str):
        """Class to store bond information

        Args:
            type1 (str): Atomtype of the first atom in the bond
            type2 (str): Atomtype of the second atom in the bond_
            bondtype (str): Type of the bond. Suported types are "harm" (harmonic) and "cons" (constraints)
            r_eq (str): Equilibrium position (A)
            k (str): Spring constant (kJ/mol)
        """
        self.atomtypes = [type1, type2]
        self.bondtype = bondtype
        self.r_eq = r_eq
        self.k = k

    @classmethod
    def parse_line(cls, line:str):
        """Function to create an object from a forcefield line

        Args:
            line (str):  Line describing a bond, following the format of CL&P

        Returns:
            Bond: An object of the Bond class
        """

        items = line.strip().split()
        assert len(items)==5, f"Found {len(items)} items in the line. Bonds are characterized by only 5 parameters"
        assert (items[2]=="harm" or items[2]=="cons"), f"Unsupported type {items[2]} was indicated for a bond. Supported types are harmonic (harm) and constraints (cons)"
        type1, type2, bondtype, r_eq, k = items
        return cls(type1, type2, bondtype, r_eq, k)

   
    def generate_line(self)->str:
        """Creates a line for the forcefield

        Returns:
            str: Forcefield line describing the bond, in the CL&P format
        """

        contents = [self.atomtypes[0], self.atomtypes[1], self.bondtype, self.r_eq, self.k]
        assert (contents[2]=="harm" or contents[2]=="cons"), f"Unsupported type {contents[2]} was indicated for a bond. Supported types are harmonic (harm) and constraints (cons)"
        line = "    ".join(contents)
        return line
    
class Angle():

    def __init__(self, type1:str, type2:str, type3:str, angletype:str, th_eq:str, k:str):
        """Class to store information about angles

        Args:
            type1 (str): Atomtype of the first atom in the angle
            type2 (str): Atomtype of the second atom in the angle
            type3 (str): Atomtype of the third atom in the angle
            angletype (str): Type of the angle. Supported types are "harm" (harmonic) and "cons" (constraints)
            th_eq (str): Equilibrium angle (deg)
            k (str): Spring constant (kJ/mol)
        """
        
        self.atomtypes = [type1, type2, type3]
        self.angletype = angletype
        self.th_eq = th_eq
        self.k = k

    @classmethod
    def parse_line(cls, line:str):
        """Generates an angle from a forcefield line

        Args:
            line (str): CL&P line describing the angle

        Returns:
            Angle: Object of the Angle class
        """

        items = line.strip().split()
        assert len(items)==6, f"Found {len(items)} items  in the line. Angles are characterized by only 6 parameters"
        assert (items[3]=="harm" or items[3]=="cons"), f"Unsupported type {items[3]} was indicated for an angle. Supported types are harmonic (harm) and constraints (cons)"
        type1, type2, type3, bondtype, th_eq, k = items
        return cls(type1, type2, type3, bondtype, th_eq, k)
    
   
    def generate_line(self)->str:
        """Generates a forcefield line from the stored information

        Returns:
            str: Forcefield line in the CL&P format
        """

        contents = [self.atomtypes[0], self.atomtypes[1], self.atomtypes[2], self.angletype, self.th_eq, self.k]
        assert (contents[3]=="harm" or contents[3]=="cons"), f"Unsupported type {contents[3]} was indicated for a bond. Supported types are harmonic (harm) and constraints (cons)"
        line = "    ".join(contents)
        return line
    

class Dihedral():

    def __init__(self, type1:str, type2:str, type3:str, type4:str, dihedraltype:str, v_params:list):
        """Class to store information about dihedrals

        Args:
            type1 (str): Atomtype of the first atom in the dihedral
            type2 (str): Atomtype of the second atom in the dihedral
            type3 (str): Atomtype of the third atom in the dihedral
            type4 (str): Atomtype of the fourth atom in the dihedral
            dihedraltype (str): Type of the dihedral. Supported types are "opls"
            v_params (list): Parameters for the dihedral. For "opls", the list contains [V1, V2, V3, V4]
        """

    
        self.atomtypes = [type1, type2, type3, type4]
        self.dihedraltype = dihedraltype
        self.params = v_params

    @classmethod
    def parse_line(cls, line:str):
        """Generates a dihedral from a forcefield line

        Args:
            line (str): CL&P line describing the dihedral

        Returns:
            Dihedral: Object of the Dihedral class
        """

        items = line.strip().split()
        assert len(items)==9, f"Found {len(items)} items  in the line. Dihedrals are characterized by only 9 parameters"
        assert items[4]=="opls", f"Unsupported type {items[4]} was indicated for a dihedral. Supported types are OPLS (opls)"
        type1, type2, type3, type4, bondtype, v1, v2, v3, v4 = items
        return cls(type1, type2, type3, type4, bondtype, [v1, v2, v3, v4])
    
   
    def generate_line(self)->str:
        """Generates a forcefield line from the stored information

        Returns:
            str: Forcefield line in the CL&P format
        """

        contents = [self.atomtypes[0], self.atomtypes[1], self.atomtypes[2], self.atomtypes[3], self.dihedraltype, self.params[0], self.params[1], self.params[2], self.params[3]]
        assert contents[4]=="opls", f"Unsupported type {contents[4]} was indicated for a dihedral. Supported types are OPLS (opls)"
        line = "    ".join(contents)
        return line


class ImproperDihedral():

    def __init__(self, type1:str, type2:str, type3:str, type4:str, dihedraltype:str, v_params:str):
        """Class to store information about improper dihedrals

        Args:
            type1 (str): Atomtype of the first atom in the dihedral
            type2 (str): Atomtype of the second atom in the dihedral
            type3 (str): Atomtype of the third atom in the dihedral
            type4 (str): Atomtype of the fourth atom in the dihedral
            dihedraltype (str): Type of the improper dihedral. Supported types are "opls"
            v_params (list): Parameters for the improper dihedral. For "opls", the list contains [V1, V2, V3, V4]
        """
        self.atomtypes = [type1, type2, type3, type4]
        self.dihedraltype = dihedraltype
        self.params = v_params

    @classmethod
    def parse_line(cls, line:str):
        """Generates an improper dihedral from a forcefield line

        Args:
            line (str): CL&P line describing the improper dihedral

        Returns:
            Dihedral: Object of the Dihedral class
        """

        items = line.strip().split()
        assert len(items)==9, f"Found {len(items)} items  in the line. Improper dihedrals are characterized by only 9 parameters"
        assert items[4]=="opls", f"Unsupported type {items[4]} was indicated for an imporper dihedral. Supported types are OPLS (opls)"
        type1, type2, type3, type4, bondtype, v1, v2, v3, v4 = items
        return cls(type1, type2, type3, type4, bondtype, [v1, v2, v3, v4])
    
   
    def generate_line(self)->str:
        """Generates a forcefield line from the stored information

        Returns:
            str: Forcefield line in the CL&P format
        """

        contents = [self.atomtypes[0], self.atomtypes[1], self.atomtypes[2], self.atomtypes[3], self.dihedraltype, self.params[0], self.params[1], self.params[2], self.params[3]]
        assert contents[4]=="opls", f"Unsupported type {contents[4]} was indicated for an improper dihedral. Supported types are OPLS (opls)"
        line = "    ".join(contents)
        return line



class Forcefield():

    def __init__(self, atoms_list:list=[], bonds_list:list=[], angles_list:list=[],
                 dihedrals_list:list=[], improper_list:list=[]):
        """Class to store a CL&P formatted forcefield

        Args:
            atoms_list (list, optional): List of Atom objects. Defaults to [].
            bonds_list (list, optional): List of Bond objects. Defaults to [].
            angles_list (list, optional): List of Angle objects. Defaults to [].
            dihedrals_list (list, optional): List of Dihedral objects. Defaults to [].
            improper_list (list, optional): List of ImproperDihedral objects. Defaults to [].
        """

        self.atoms = atoms_list
        self.bonds = bonds_list
        self.angles = angles_list
        self.dihedrals = dihedrals_list
        self.improper = improper_list


    @classmethod
    def from_file(cls, forcefield_file:str):
        """Generates a forcefield from a file

        Args:
            forcefield_file (str): Path to a file containing a CL&P forcefield

        Returns:
            Forcefield: Object of the Forcefield class
        """

        with open(forcefield_file, "r") as file:

            lines = file.readlines()

        SECTIONS = {"ATOMS":[], "BONDS":[], "ANGLES":[], "DIHEDRALS":[], "IMPROPER":[]}
        CLASSES = {"ATOMS":Atom, "BONDS":Bond, "ANGLES":Angle, "DIHEDRALS":Dihedral, "IMPROPER":ImproperDihedral}
        CURRENT_SECTION = None
        for line in lines:

            line = line.strip()

            if not line.startswith("#") and line != "":

                if line in SECTIONS.keys():
                    CURRENT_SECTION = line
                else:
                    SECTIONS[CURRENT_SECTION].append(CLASSES[CURRENT_SECTION].parse_line(line))
        
        return cls(SECTIONS["ATOMS"], SECTIONS["BONDS"], SECTIONS["ANGLES"], SECTIONS["DIHEDRALS"], SECTIONS["IMPROPER"])

                
    def get_sub_forcefield(self, atom_names:list):
        """Generates a subset of the forcefield that only contains information about
           atoms of interest

        Args:
            atom_names (list): List with the names of the atoms which will appear in the sub-forcefield

        Returns:
            Forcefield: A new object of the Forcefield class, which contains information only about the 
                        given atom_names (and the bonds, angles,... etc. that involve only the selected atoms).
        """
        
        atoms_list = []
        bonds_list = []
        angles_list = []
        dihedrals_list = []
        improper_list = []

        ## We select only the atoms named in the selection
        for atom in self.atoms:
            if atom.name in atom_names:
                atoms_list.append(atom)

        ## Check if we missed any atoms
        _selected_names = [a.name for a in atoms_list]
        for atomname in atom_names:
            assert atomname in _selected_names, f"Atom with name {atomname} is not in the forcefield"

        ## Now we get the corresponding atomtypes
        atomtypes = set([a.atomtype for a in atoms_list])

        ## We now select only bonds angles and dihedrals from the original forcefield
        ## if and only if every atomtype they involve is in the selection

        for bond in self.bonds:
            if set(bond.atomtypes).issubset(atomtypes):
                bonds_list.append(bond)

        for angle in self.angles:
            if set(angle.atomtypes).issubset(atomtypes):
                angles_list.append(angle)
        
        for dihedral in self.dihedrals:
            if set(dihedral.atomtypes).issubset(atomtypes):
                dihedrals_list.append(dihedral)

        for improper in self.improper:
            if set(improper.atomtypes).issubset(atomtypes):
                improper_list.append(improper)

        return Forcefield(atoms_list, bonds_list, angles_list, dihedrals_list, improper_list)
        
    def write_file(self, filename:str, comment:str="## Written with metaltool"):
        """Writes the current forcefield information

        Args:
            filename (str): Name (or path) of the new forcefield file
            comment (str, optional): Comment line in the forcefield. Defaults to "## Written with metaltool".
        """
        with open(filename, "w") as file:

            file.write(comment)

            if len(self.atoms) != 0:
                file.write("\nATOMS\n")
                file.write("#     type   m/u     q/e    pot    *parameters    polarizability/A-3\n")
                for atom in self.atoms:
                    file.write(atom.generate_line())
                    file.write("\n")

            if len(self.bonds) != 0:
                file.write("\nBONDS\n")
                file.write("\n# i   j    pot    re/A    kr/kJmol-1\n")
                for bond in self.bonds:
                    file.write(bond.generate_line())
                    file.write("\n")

            if len(self.angles) != 0:
                file.write("\nANGLES\n")
                file.write("# i   j   k    pot    th/deg  ka/kjmol-1\n")
                for angle in self.angles:
                    file.write(angle.generate_line())
                    file.write("\n")

            if len(self.dihedrals) != 0:
                file.write("\nDIHEDRALS\n")
                file.write("# i   j   k   l    pot     v1        v2        v3        v4\n")
                for dihedral in self.dihedrals:
                    file.write(dihedral.generate_line())
                    file.write("\n")

            if len(self.improper) != 0:
                file.write("\nIMPROPER\n")
                file.write("# i   j   k   l    pot     v1        v2        v3        v4\n")
                for improper in self.improper:
                    file.write(improper.generate_line())
                    file.write("\n")

def get_atomnames(file:str)->list:
    """Given a structure file, retrieves the names of the atoms present in that file

    Args:
        file (str): Path to the structure file

    Raises:
        ValueError: If a not supported structure format is provided. Supported formats are zmat and xyz

    Returns:
        list: List containing the names of the atoms in the provided file
    """

    extension = os.path.splitext(file)[-1]
    atomnames = []

    if extension==".zmat":

        with open(file, "r") as f:
            lines = f.readlines()
        
        for line in lines[2:]:
            items = line.strip().split()

            ## We break at the first blank line
            try:
                atomnames.append(items[1])
            except:
                break


    elif extension==".xyz":
        
        with open(file, "r") as f:
            lines = f.readlines()

        for line in lines[2:]:
            atomnames.append(line.strip().split()[0])


    else:
        raise ValueError(f"Files with extension {extension} are not supported. Supported extensions are .zmat and .xyz")
    
    return atomnames


def main(forcefield: Forcefield, files: list):
    """
    Given a forcefield and a list of structure files, generates input files for fftool
    where no atoms have repeating names and only the relevant parameters will appear in
    the resulting forcefield

    Args:
        forcefield (Forcefield): Object of the Forcefield class
        file (str): List of paths to .zmat or .xyz containing structure files
    """

    original_atomnames = []
    for file in files:
        names = get_atomnames(file)
        original_atomnames += names


    sub_ff = forcefield.get_sub_forcefield(set(original_atomnames))

    ## Now we need to duplicate some atoms in the new forcefield, since repeated atoms will get a new name, 
    ## and thus new entries will be required in the forcefield (even though the new different names will have the same atomtype)

    repetitions = {key:number for key, number in Counter(original_atomnames).items() if number > 1}

    for atomname, count in repetitions.items():

        atom = [a for a in sub_ff.atoms if a.name==atomname][0]

        ## We will add copies of the repeated atoms to the forcefield. For now, they will be identical to one another
        ## Later the names will be changed
        for _ in range(count-1):

            new_atom = Atom(atom.name, atom.atomtype, atom.mass, atom.charge, atom.pottype, atom._potential_params, atom.polarizability)
            sub_ff.atoms.append(new_atom)

    ## Now we will open the old zmat/xyz files to create new ones, with the atoms enumerated in ascending order,
    ## this way there won't be repeated atoms.
    
    atom_id = 1
    for file in files:
        extension = os.path.splitext(file)[-1]

        if extension==".zmat":

            new_lines = []

            with open(file, "r") as f:
                original_lines = f.readlines()

            new_lines.append(original_lines[0])
            new_lines.append(original_lines[1])

            _line_index = 2
            for line in original_lines[2:]:
                items = line.split()
                try:
                    ## Change the name of the atom in the file
                    old_name = items[1]
                    new_name = items[1][:2].rstrip("0123456789") + str(atom_id)
                    items[1] = new_name
                    new_lines.append("    ".join(items)+"\n")
                    ## Change the name of the atom in the forcefield
                    [a for a in sub_ff.atoms if a.name==old_name][0].name = new_name
                    ## Increase the atom id for the newt line, as well as the line index
                    atom_id += 1
                    _line_index += 1
                except:
                    break

            ## Add the rest of the lines to the new file
            
            new_lines += original_lines[_line_index:]
            new_lines[-1] = "new_forcefield.ff"

            ## Save the new file

            with open(f"new_{file}", "w") as f:
                f.writelines(new_lines)


        elif extension==".xyz":

            new_lines = []
            with open(file, "r") as f:
                original_lines = f.readlines()

            new_lines.append(original_lines[0])
            molecule, _ = original_lines[1].split()
            new_lines.append(f"{molecule} new_forcefield.ff\n")

            for line in original_lines[2:]:
                items = line.split()
                old_name = items[0]
                new_name = items[0][:2].rstrip("0123456789") + str(atom_id)
                items[0] = new_name
                new_lines.append("      ".join(items)+"\n")
                ## Change the name of the atom in the forcefield
                print([a for a in sub_ff.atoms if a.name==old_name][0].name)
                [a for a in sub_ff.atoms if a.name==old_name][0].name = new_name
                ## Increase the atom id for the newt line
                atom_id += 1

            ## Save the new file

            with open(f"new_{file}", "w") as f:
                f.writelines(new_lines)
        
        else:
            raise ValueError(f"File extension {extension} is not supported. Supported extensions are zmat and xyz")
    ## Save the forcefield
                
    sub_ff.write_file("new_forcefield.ff")

        
if __name__ == '__main__':

    forcefield_file = sys.argv[1]
    structure_files = sys.argv[2:]
    initial_ff = Forcefield.from_file(forcefield_file)
    main(initial_ff, structure_files)



