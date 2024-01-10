import sys
from collections import Counter
import os


class Atom():

    def __init__(self, name, type, mass, charge,
                 potential_type, potential_params, polarizability="0.0"):

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
    def parse_line(cls, line):

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

    def generate_line(self):

        if self.pottype == "lj":

            contents = [self.name, self.atomtype, self.mass, self.charge, self.pottype, self.sigma, self.epsilon, self.polarizability]

        elif self.pottype == "ft":

            contents = [self.name, self.atomtype, self.mass, self.charge, self.pottype, self.B, self.A, self.C6, self.C8, self.B, self.B, self.polarizability]

        else:
            raise ValueError(f"Atom {self.name} of type {self.atomtype} has been assigned a {self.pottype} potential, which is not supported. Supported types are Lennard-Jones (lj) and Fumi-Tosi (ft)")

        line = "    ".join(contents)

        return line
    

class Bond():

    def __init__(self, type1, type2, bondtype, r_eq, k):

        self.atomtypes = [type1, type2]
        self.bondtype = bondtype
        self.r_eq = r_eq
        self.k = k

    @classmethod
    def parse_line(cls, line):

        items = line.strip().split()
        assert len(items)==5, f"Found {len(items)} items in the line. Bonds are characterized by only 5 parameters"
        assert (items[2]=="harm" or items[2]=="cons"), f"Unsupported type {items[2]} was indicated for a bond. Supported types are harmonic (harm) and constraints (cons)"
        type1, type2, bondtype, r_eq, k = items
        return cls(type1, type2, bondtype, r_eq, k)

   
    def generate_line(self):

        contents = [self.atomtypes[0], self.atomtypes[1], self.bondtype, self.r_eq, self.k]
        assert (contents[2]=="harm" or contents[2]=="cons"), f"Unsupported type {contents[2]} was indicated for a bond. Supported types are harmonic (harm) and constraints (cons)"
        line = "    ".join(contents)
        return line
    
class Angle():

    def __init__(self, type1, type2, type3, angletype, th_eq, k):
        
        self.atomtypes = [type1, type2, type3]
        self.angletype = angletype
        self.th_eq = th_eq
        self.k = k

    @classmethod
    def parse_line(cls, line):

        items = line.strip().split()
        assert len(items)==6, f"Found {len(items)} items  in the line. Angles are characterized by only 6 parameters"
        assert (items[3]=="harm" or items[3]=="cons"), f"Unsupported type {items[3]} was indicated for an angle. Supported types are harmonic (harm) and constraints (cons)"
        type1, type2, type3, bondtype, th_eq, k = items
        return cls(type1, type2, type3, bondtype, th_eq, k)
    
   
    def generate_line(self):

        contents = [self.atomtypes[0], self.atomtypes[1], self.atomtypes[2], self.angletype, self.th_eq, self.k]
        assert (contents[3]=="harm" or contents[3]=="cons"), f"Unsupported type {contents[3]} was indicated for a bond. Supported types are harmonic (harm) and constraints (cons)"
        line = "    ".join(contents)
        return line
    

class Dihedral():

    def __init__(self, type1, type2, type3, type4, dihedraltype, v_params):
    
        self.atomtypes = [type1, type2, type3, type4]
        self.dihedraltype = dihedraltype
        self.params = v_params

    @classmethod
    def parse_line(cls, line):

        items = line.strip().split()
        assert len(items)==9, f"Found {len(items)} items  in the line. Dihedrals are characterized by only 9 parameters"
        assert items[4]=="opls", f"Unsupported type {items[4]} was indicated for a dihedral. Supported types are OPLS (opls)"
        type1, type2, type3, type4, bondtype, v1, v2, v3, v4 = items
        return cls(type1, type2, type3, type4, bondtype, [v1, v2, v3, v4])
    
   
    def generate_line(self):

        contents = [self.atomtypes[0], self.atomtypes[1], self.atomtypes[2], self.atomtypes[3], self.dihedraltype, self.params[0], self.params[1], self.params[2], self.params[3]]
        assert contents[4]=="opls", f"Unsupported type {contents[4]} was indicated for a dihedral. Supported types are OPLS (opls)"
        line = "    ".join(contents)
        return line


class ImproperDihedral():

    def __init__(self, type1, type2, type3, type4, dihedraltype, v_params):
    
        self.atomtypes = [type1, type2, type3, type4]
        self.dihedraltype = dihedraltype
        self.params = v_params

    @classmethod
    def parse_line(cls, line):

        items = line.strip().split()
        assert len(items)==9, f"Found {len(items)} items  in the line. Improper dihedrals are characterized by only 9 parameters"
        assert items[4]=="opls", f"Unsupported type {items[4]} was indicated for an imporper dihedral. Supported types are OPLS (opls)"
        type1, type2, type3, type4, bondtype, v1, v2, v3, v4 = items
        return cls(type1, type2, type3, type4, bondtype, [v1, v2, v3, v4])
    
   
    def generate_line(self):

        contents = [self.atomtypes[0], self.atomtypes[1], self.atomtypes[2], self.atomtypes[3], self.dihedraltype, self.params[0], self.params[1], self.params[2], self.params[3]]
        assert contents[4]=="opls", f"Unsupported type {contents[4]} was indicated for an improper dihedral. Supported types are OPLS (opls)"
        line = "    ".join(contents)
        return line



class Forcefield():

    def __init__(self, atoms_list=[], bonds_list=[], angles_list=[],
                 dihedrals_list=[], improper_list=[]):

        self.atoms = atoms_list
        self.bonds = bonds_list
        self.angles = angles_list
        self.dihedrals = dihedrals_list
        self.improper = improper_list


    @classmethod
    def from_file(cls, forcefield_file):

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

                
    def get_sub_forcefield(self, atom_names):
        
        atoms_list = []
        bonds_list = []
        angles_list = []
        dihedrals_list = []
        improper_list = []

        ## We select only the atoms named in the selection
        for atom in self.atoms:
            if atom.name in atom_names:
                atoms_list.append(atom)

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
        
    def write_file(self, filename, comment="## Written with metaltool"):

        with open(filename, "w") as file:

            file.write(comment)
            file.write("\n")

            if len(self.atoms) != 0:
                file.write("ATOMS\n")
                for atom in self.atoms:
                    file.write(atom.generate_line())
                    file.write("\n")

            if len(self.bonds) != 0:
                file.write("BONDS\n")
                for bond in self.bonds:
                    file.write(bond.generate_line())
                    file.write("\n")

            if len(self.angles) != 0:
                file.write("ANGLES\n")
                for angle in self.angles:
                    file.write(angle.generate_line())
                    file.write("\n")

            if len(self.dihedrals) != 0:
                file.write("DIHEDRALS\n")
                for dihedral in self.dihedrals:
                    file.write(dihedral.generate_line())
                    file.write("\n")

            if len(self.improper) != 0:
                file.write("IMPROPER\n")
                for improper in self.improper:
                    file.write(improper.generate_line())
                    file.write("\n")

def get_atomnames(file):

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
                ## Increase the atom id for the newt line, as well as the line index
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



