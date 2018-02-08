
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
    contact_single_regex = re.compile("^([A-Za-z]+[']?)([0-9]+)$")
    contact_pair_regex_1 = re.compile("^([A-Za-z]+[']?)([0-9]+)-([0-9]+)$")
    contact_pair_regex_2 = re.compile("^([A-Za-z]+[']?)([0-9]+)-([A-Za-z]+[']?)([0-9]+)$")

    def __init__(self, filename):
        super().__init__(filename)
        self.electrodes = {}
        for i, name in enumerate(self.names):
            match = self.contact_single_regex.match(name)
            if match is None:
                raise ValueError("Unexpected contact name %s" % name)

            elec_name, _ = match.groups()
            if elec_name not in self.electrodes:
                self.electrodes[elec_name] = []
            self.electrodes[elec_name].append(i)

    def get_elec(self, name):
        match = self.contact_single_regex.match(name)
        if match is None:
            return None

        return match.groups()[0]

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
        """

        match = self.contact_single_regex.match(name)
        if match is not None:
            return self.name_to_xyz[name]

        match = self.contact_pair_regex_1.match(name)
        if match is not None:
            assert abs(int(match.group(2)) - int(match.group(3))) == 1
            contact1 = match.group(1) + match.group(2)
            contact2 = match.group(1) + match.group(3)
            return (self.name_to_xyz[contact1] + self.name_to_xyz[contact2])/2.

        match = self.contact_pair_regex_2.match(name)
        if match is not None:
            assert match.group(1) == match.group(3)
            assert abs(int(match.group(2)) - int(match.group(4))) == 1
            contact1 = match.group(1) + match.group(2)
            contact2 = match.group(3) + match.group(4)
            return (self.name_to_xyz[contact1] + self.name_to_xyz[contact2])/2.

        raise ValueError("Given name '%s' does not follow any expected pattern." % name)
