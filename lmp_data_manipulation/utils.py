class f_list:
    """
    WIP

    List of floats

    """

    def __init__(self, elements=[]):
        self.elements = elements

    def __len__(self):
        return len(self.elements)

    def __str__(self):
        return '[{0}]'.format(', '.join(str(el) for el in self.elements))

    def __getitem__(self, idx):
        return self.elements[idx]


    def ave(self):
        if not self.elements:
            return None
        return sum(self.elements) / len(self.elements)


def ae(a, b, epsilon=0.01):
    return abs(a-b) < epsilon*(abs(a) + abs(b))
