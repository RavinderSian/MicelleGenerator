import numpy as np

class Micelle():

    def __init__(self):
        self.user_input = input("[filename]  [number of atoms]")  # take in .gro file and number of atoms
        self.molecules = int(input("how many molecules?"))  # take in number of molecules and set as variable
        self.head_atom = int(input("What is the head atom number?"))  # takes in the head atom
        self.no_of_atoms = int(self.user_input.split()[1])  # Makes number of atoms equal to index 1 of userinput which is the number of atoms
        file_name = self.user_input.split()[0]  # file name input
        self.coords = [coord.rstrip("\n") for coord in open(file_name)]  # removes “/n” literal, which symbolises a new line from the readout of the .gro file within the script

    def move_to_surface(self):  # moves one single lipopeptide to micelle surface
        new_list = []  # creates empty list for storage of translated coordinates
        head_coord = self.coords[5].split()[3:]  # turns the head group coordinates into a variable
        second_line = self.molecules * self.no_of_atoms  # obtains total number of atoms as a variable
        radius = float(input("enter the micelle radius in nanometers"))  # takes in radius input from user
        radii = str(radius)  # converts radius to string for printing into text file
        opening_lines = "micelle of radius {:>1} nm \n {:>1}".format(radius, second_line)  # generates opening lines of .gro output file
        self.new_name = "micelle" + radii + "nm" "molecules" + str(self.molecules) + ".gro"  # defines file name of output in .gro format
        print(opening_lines, file=open(self.new_name, "a"))  # creates output file and prints title and number of atoms to output file
        line_number = 2  # sets line index to 2 so first 2 lines are ignored

        while line_number < self.no_of_atoms + 2:  # while loop performs translation of single lipopeptide line by line
            elements = self.coords[line_number].split()  # turns line of coordinates for index line_number into a list
            head_list = [float(head_coord[0]), float(head_coord[1]), float(head_coord[2])]  # creates a list of headgroup coordinates
            head_arr = np.array(head_list)  # turns head_list to array
            x_coordinate = float(elements[3])
            y_coordinate = float(elements[4])  # turns string coordinates into numbers
            z_coordinate = float(elements[5]) + radius  # translates lipopeptide by radius input along z axis only
            coord_list = [x_coordinate, y_coordinate, z_coordinate]  # turns coordinates of index line_number into a list of numbers
            coord_arr = np.array(coord_list)  # list of coordinates in number form is made into an array
            coord_surface = coord_arr - head_arr  # translates current coordinates by head group coordinates
            new_list.append(coord_surface)  # array of translated coordinates is appended to an empty list
            self.array = new_list  # the list is made into a variable so it can be made into an array and rotated
            line_number += 1  # adds one to line_number counter for each time the while loop is run

    def rotation(self):  # copies and rotates lipopeptide at micelle surface
        molecules_count = 0  # sets molecule counter to 0
        self.rot_list = []  # creates empty list for rotated coordinates
        while self.molecules > molecules_count:
            molecules_count += 1  # adds 1 to molecule counter
            rand_numb_theta = np.random.uniform(low=-1)  # generates random theta value
            random_radian_theta = np.arccos(rand_numb_theta)  # returns inverse cos of theta generated in radians
            random_theta = np.rad2deg(random_radian_theta)  # converts random theta value from radians to degrees
            phi = np.random.uniform(low=-1) * 180  # generates value for phi
            line_number = 0  # sets line number counter to 0, acts as loop counter
            while line_number < self.no_of_atoms:  # loop ends when number of lines is equal to number of atoms
                array_for_rotation = np.eye(3)  # generates a 3x3 array with 1 diagonal and zeros in all other poisitions
                array_for_rotation_by_x_axis = np.copy(array_for_rotation)  # creates basis for rotation matrix for rotation around x axis
                array_for_rotation_by_z_axis = np.copy(array_for_rotation)  # creates basis for rotation matrix for rotation around z axis
                self.rot = np.array(self.array)[line_number].reshape(3, 1)
                line_number += 1  # adds 1 to loop counter

                array_for_rotation_by_x_axis[1][1] = np.cos(np.deg2rad(random_theta))
                array_for_rotation_by_x_axis[1][2] = -np.sin(np.deg2rad(random_theta))  # rotation matrix for rotation around x axis
                array_for_rotation_by_x_axis[2][1] = np.sin(np.deg2rad(random_theta))
                array_for_rotation_by_x_axis[2][2] = np.cos(np.deg2rad(random_theta))
                rotbyxmat = np.asmatrix(array_for_rotation_by_x_axis).round(3)  # turns array into matrix

                array_for_rotation_by_z_axis[0][0] = np.cos(np.deg2rad(phi))
                array_for_rotation_by_z_axis[0][1] = -np.sin(np.deg2rad(phi))  # rotation matrix for rotation by z axis
                array_for_rotation_by_z_axis[1][0] = np.sin(np.deg2rad(phi))
                array_for_rotation_by_z_axis[1][1] = np.cos(np.deg2rad(phi))
                rotbyzmat = np.asmatrix(array_for_rotation_by_z_axis).round(3)

                longitude = np.matmul(rotbyxmat, self.rot)  # coordinates rotated around x axis
                latitude = np.matmul(rotbyzmat, longitude).round(3)  # coordinates rotated around z and x axis
                final_latitude = latitude.reshape(1, 3)
                self.rot_list.append(final_latitude)

    def formatting(self):
        atom_number = 1  # sets atom number to 1
        line_counter = 0  # loop counter
        lineno_list = 0  # used to end loop
        residue_number = 1  # sets residue number to 1
        line_number = 0  # sets line number to 0
        residue = self.coords[2:][0][5:9]  # defines current residue
        next_residue = self.coords[2:][1][5:9]  # defines residue of next line

        while line_counter < self.molecules * self.no_of_atoms:
            if lineno_list < self.no_of_atoms:
                x = self.rot_list[line_counter][0][0]
                y = self.rot_list[line_counter][0][1]  # uses a while loop and print output file line by line
                z = self.rot_list[line_counter][0][2]
                rounded_x_coordinate = "%.3f" % x
                rounded_y_coordinate = "%.3f" % y  # keeps to 3 decimal places
                rounded_z_coordinate = "%.3f" % z
                atom_type = self.coords[2:][lineno_list].split()[1]
                self.coordinates = "{:>3}{:>5}{:>8}{:>8}{:>8}".format(atom_type, atom_number, rounded_x_coordinate, rounded_y_coordinate, rounded_z_coordinate)
                # defines new line for output file
                atom_number += 1
                line_counter += 1
                lineno_list += 1
            if lineno_list == self.no_of_atoms:
                lineno_list = 0
            if line_number < self.no_of_atoms - 1 and residue == next_residue:
                residue = self.coords[2:][line_number][5:9]
                next_residue = self.coords[2:][line_number + 1][5:9]
                final_line = "{:>5}{:>4}{:>35}".format(residue_number, residue, self.coordinates)  # keeps same residue for next line
                print(final_line, file=open(self.new_name, "a"))
                line_number += 1
            elif line_number == self.no_of_atoms - 1:
                residue = self.coords[2:][self.no_of_atoms - 1][5:9]
                next_residue = self.coords[2:][0][5:9]

            while line_counter < self.molecules * self.no_of_atoms:
                if lineno_list < self.no_of_atoms:
                    x = self.rot_list[line_counter][0][0]
                    y = self.rot_list[line_counter][0][1]  # uses a while loop and print output file line by line
                    z = self.rot_list[line_counter][0][2]
                    rounded_x_coordinate = "%.3f" % x
                    rounded_y_coordinate = "%.3f" % y  # keeps to 3 decimal places
                    rounded_z_coordinate = "%.3f" % z
                    atom_type = self.coords[2:][lineno_list].split()[1]
                    self.coordinates = "{:>3}{:>5}{:>8}{:>8}{:>8}".format(atom_type, atom_number, rounded_x_coordinate, rounded_y_coordinate, rounded_z_coordinate)  # defines new line for output file
                    atom_number += 1
                    line_counter += 1
                    lineno_list += 1
                if lineno_list == self.no_of_atoms:
                    lineno_list = 0
                if line_number < self.no_of_atoms - 1 and residue == next_residue:
                    residue = self.coords[2:][line_number][5:9]
                    next_residue = self.coords[2:][line_number + 1][5:9]
                    final_line = "{:>5}{:>4}{:>35}".format(residue_number, residue, self.coordinates)  # keeps same residue for next line
                    print(final_line, file=open(self.new_name, "a"))
                    line_number += 1
                elif line_number == self.no_of_atoms - 1:
                    residue = self.coords[2:][self.no_of_atoms - 1][5:9]
                    next_residue = self.coords[2:][0][5:9]
                    final_line = "{:>5}{:>4}{:>35}".format(residue_number, residue, self.coordinates)  # resets loop when lipopeptide number of atoms is exceeded
                    print(final_line, file=open(self.new_name, "a"))
                    line_number = 0
                elif line_number < self.no_of_atoms - 1 and residue != next_residue:  # adds 1 to residue number if next line has different residue to current line
                    residue_number += 1
                    residue = self.coords[2:][line_number][5:9]
                    next_residue = self.coords[2:][line_number + 1][5:9]
                    final_line = "{:>5}{:>4}{:>35}".format(residue_number, residue, self.coordinates)
                    print(final_line, file=open(self.new_name, "a"))
                    line_number += 1
        print("1.0 1.0 1.0", file=open(self.new_name, "a"))


micelle: Micelle = Micelle()
micelle.move_to_surface()
micelle.rotation()
micelle.formatting()
