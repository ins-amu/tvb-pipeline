
import io
import re

import numpy as np


class NamedPoints():
    def __init__(self, fl):
        data = np.genfromtxt(fl, dtype=None)
        self.xyz = np.array([[l[1], l[2], l[3]] for l in data])
        self.names = [l[0].decode('ascii') for l in data]
        self.name_to_xyz = dict(zip(self.names, self.xyz))


class Contacts(NamedPoints):
    """
    Load the contacts names and position from file.
    Assumptions:
    - Numbering of every electrode starts from 1 and is contiguous: 1, 2, 3, 4, ...
    - Number of the contact is not prepended by zero.
    - Name of the electrode can contain letters (upper- and lowercase), numbers, and "'".

    >>> contacts = Contacts(io.BytesIO("A1 0 0 1\\nA2 0 0 2".encode()))
    >>> list(contacts.electrodes.keys())
    ['A']

    >>> contacts = Contacts(io.BytesIO("A11 0 0 1\\nA12 0 0 2\\nA21 1 0 1\\nA22 1 0 2".encode()))
    >>> list(contacts.electrodes.keys())
    ['A1', 'A2']
    """
    contact_pair_regex = re.compile("^([A-Za-z0-9']+)-([A-Za-z0-9']+)$")

    def __init__(self, filename):
        super().__init__(filename)
        self.electrodes = {}
        self.name_to_elec = {}

        i = 0
        while i < len(self.names):
            if self.names[i][-1] != "1":
                raise ValueError("Numbering of contacts on electrodes should start at 1. (Line %d, '%s')" % (i, name))
            elec_name = self.names[i][:-1]
            elec_num = 0
            self.electrodes[elec_name] = []

            while i < len(self.names):
                if self.names[i] != elec_name + str(elec_num + 1):
                    break

                self.electrodes[elec_name].append(i)
                self.name_to_elec[self.names[i]] = elec_name
                elec_num = elec_num + 1
                i = i +1


    def get_elec(self, name):
        return self.name_to_elec.get(name, None)


    def get_coords(self, name):
        """Get the coordinates of a specified contact or contact pair. Allowed formats are:

        A1            : Single contact.
        A1-2 or A1-A2 : Contact pair. The indices must be adjacent.


        Examples:

        >>> np.set_printoptions(formatter={'float': lambda x: "{0:0.1f}".format(x)})

        >>> contacts = Contacts(io.BytesIO("A1 0.0 0.0 1.0\\nA2 0.0 0.0 2.0".encode()))
        >>> contacts.get_coords("A1")
        array([0.0, 0.0, 1.0])

        >>> contacts.get_coords("A1-2")
        array([0.0, 0.0, 1.5])

        >>> contacts.get_coords("A2-A1")
        array([0.0, 0.0, 1.5])

        >>> contacts.get_coords("Invalid")
        Traceback (most recent call last):
        ...
        ValueError: Cannot get coordinates for 'Invalid'.
        """


        # Solo contact
        if name in self.names:
            return self.name_to_xyz[name]

        # Contact pair
        match = self.contact_pair_regex.match(name)
        if match is not None:
            a1, a2 = match.group(1), match.group(2)

            # Model A1-A2
            if a1 in self.names and a2 in self.names:
                return (self.name_to_xyz[a1] + self.name_to_xyz[a2])/2.

            # Model A1-2
            elif a1 in self.names:
                a2 = self.name_to_elec[a1] + a2
                if a2 in self.names:
                    return (self.name_to_xyz[a1] + self.name_to_xyz[a2])/2.

        raise ValueError("Cannot get coordinates for '%s'." % name)
