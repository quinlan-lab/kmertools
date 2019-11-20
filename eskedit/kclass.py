from collections import defaultdict


class Variant:
    def __init__(self, *args, **kwargs):
        if len(kwargs) == 0:
            self.REF = args[0]
            self.ALT = args[1]
            self.POS = args[2]
            self.CHROM = args[3]
        else:
            variant = kwargs.get('variant')
            try:
                fields = kwargs.get('fields')
            except KeyError:
                fields = ['vep']
            self.REF = variant.REF
            self.ALT = variant.ALT
            self.POS = variant.POS
            self.CHROM = variant.CHROM
            self.info = defaultdict(str)
            for f in fields:
                self.info[f] = str(variant.INFO.get(f))
            self.AC = variant.INFO.get('AC')
        self.INDEX = self.CHROM + ' ' + str(self.POS)

    def __str__(self):
        # return "POS: " + str(self.POS) + "\tREF: " + str(self.REF) + "\tALT: " + str(self.ALT) + '\n'
        return str(self.CHROM) + "\t" + str(self.POS) + "\t" + str(self.REF) + "\t" + str(self.ALT) + "\t" + str(
            self.AC) + '\n'

    def __repr__(self):
        return "CHROM: " + str(self.CHROM) + "\t" + "POS: " + str(self.POS) + "\tREF: " + str(
            self.REF) + "\tALT: " + str(self.ALT) + "\t" + str(self.AC) + '\n'

    def __eq__(self, other):
        if not isinstance(other, Variant):
            return False
        return other.CHROM == self.CHROM and other.POS == self.POS and other.REF == self.REF and other.ALT == self.ALT

    def __ne__(self, other):
        if not isinstance(other, Variant):
            return True
        return not (
                other.CHROM == self.CHROM and other.POS == self.POS and other.REF == self.REF and other.ALT == self.ALT)

    def __hash__(self):
        return hash(str(self))

    def __lt__(self, other):
        if other.CHROM != self.CHROM:
            return self.CHROM < other.CHROM
        else:
            return self.POS < other.POS

    def __gt__(self, other):
        if other.CHROM != self.CHROM:
            return self.CHROM > other.CHROM
        else:
            return self.POS > other.POS

    def csv_str(self):
        return str(self.CHROM) + "\t" + str(self.POS) + "\t" + str(self.REF) + "\t" + str(self.ALT) + "\t" + str(
            self.AC) + '\n'

    def print_variant(self, fields=None):
        if fields is None:
            fields = []
        pfields = [f for f in fields if f in self.info.keys()]
        output = str(self)[:-1] + '\t'
        for i, f in enumerate(pfields):
            output += str(self.info[f])
            if i != (len(pfields) - 1):
                output += '\t'
        return output + '\n'


class VCFRegion:
    def __init__(self, chrom, start, stop):
        self.region = [chrom, start, stop]

    # def add(self, chrom, start, stop):
    #     if chrom in self.regions.keys():
    #         self.regions[chrom][1] = stop

    def size(self):
        return self.region[2] - self.region[1]

    def __str__(self):
        return str(self.region[0]) + ':' + str(self.region[1]) + '-' + str(self.region[2])

    def __repr__(self):
        return str(self)

    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other):
        return str(other) == str(self)
