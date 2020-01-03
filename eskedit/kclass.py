from collections import defaultdict

from eskedit.constants import get_grch38_chroms


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
                self.qual = kwargs.get('qual')
            except KeyError:
                self.qual = False
            try:
                fields = kwargs.get('fields')
            except KeyError:
                fields = ['vep']
            self.REF = variant.REF
            self.ALT = variant.ALT
            self.POS = variant.POS
            self.CHROM = variant.CHROM
            if self.qual:
                self.QUAL = variant.QUAL
            self.info = defaultdict(str)
            for f in fields:
                self.info[f] = str(variant.INFO.get(f))
            self.AC = variant.INFO.get('AC')
        self.INDEX = self.CHROM + ' ' + str(self.POS)

    def __str__(self):
        # return "POS: " + str(self.POS) + "\tREF: " + str(self.REF) + "\tALT: " + str(self.ALT) + '\n'
        if self.qual:
            return str(self.CHROM) + "\t" + str(self.POS) + "\t" + str(self.REF) + "\t" + str(self.ALT) + "\t" + str(
                self.AC) + '\t' + str(self.QUAL) + '\n'
        return str(self.CHROM) + "\t" + str(self.POS) + "\t" + str(self.REF) + "\t" + str(self.ALT) + "\t" + str(
            self.AC) + '\n'

    def __repr__(self):
        # return "CHROM: " + str(self.CHROM) + "\t" + "POS: " + str(self.POS) + "\tREF: " + str(
        #     self.REF) + "\tALT: " + str(self.ALT) + "\t" + str(self.AC) + '\n'
        return str(self)

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


class Kmer:
    def __init__(self, sequence):
        from eskedit import get_complementary_sequence
        self.sequence = sequence.upper()
        self.complement = get_complementary_sequence(self.sequence)
        if self.complement < self.sequence:
            temp = self.sequence
            self.sequence = self.complement
            self.complement = temp

    def __eq__(self, other):
        return other.sequence == self.sequence or other.sequence == self.complement

    def __ne__(self, other):
        return other.sequence != self.sequence and other.sequence != self.complement

    def __str__(self):
        return self.sequence

    def __repr__(self):
        return self.sequence

    def __hash__(self):
        return hash(self.sequence)

    def __lt__(self, other):
        return other.sequence > self.sequence

    def __gt__(self, other):
        return other.sequence < self.sequence

    def __getitem__(self, item):
        if not isinstance(item, (int, float)):
            return ""
        return self.sequence[item]


class VCFRegion:
    def __init__(self, chrom, start, stop):
        # self.region = [chrom, start, stop]
        self.chrom = chrom
        try:
            self.start = int(start)
            self.stop = int(stop)
        except ValueError:
            raise ValueError('Start and Stop positions must be integers. Read %s as \'start\' and %s as \'stop\'' % (
                str(start), str(stop)))

    # def add(self, chrom, start, stop):
    #     if chrom in self.regions.keys():
    #         self.regions[chrom][1] = stop
    def is_complete(self):
        if None in [self.chrom, self.start, self.stop]:
            return False
        elif int(self.stop) - int(self.start) < 0:
            return False
        else:
            return True

    def size(self):
        return self.start - self.stop

    def __str__(self):
        return str(self.chrom) + ':' + str(self.start) + '-' + str(self.stop)

    def __repr__(self):
        return str(self)

    def __hash__(self):
        return hash(str(self.chrom) + str(self.start) + str(self.stop))

    def __eq__(self, other):
        return str(other) == str(self)

    def __lt__(self, other):
        if self.chrom == other.chrom:
            return self.start < other.start
        else:
            return self.chrom < other.chrom

    # def __getstate__(self):
    #     return self.__hash__
    #
    # def __setstate__(self, state):
    #     self.__dict__ = state


class RegionContainer:
    def __init__(self):
        # Store as a dictionary indexed by chromosome, mapping to regions
        self.regions = defaultdict(list)

    def add_distinct_region(self, region):
        current_regions = self.regions[region.chrom]
        current_regions.append(region)
        self.regions[region.chrom] = current_regions

    def add_region(self, region):
        current_regions = set(self.regions[region.chrom])
        modified = False
        for sreg in current_regions:
            if region.start < sreg.start < region.stop:
                if region.stop > sreg.stop:
                    sreg.start = region.start
                    sreg.stop = region.stop
                    modified = True
                elif region.stop <= sreg.stop:
                    sreg.start = region.start
                    modified = True
            elif (sreg.start < region.start < sreg.stop) and (region.stop <= sreg.stop):
                modified = True
            elif sreg.start < region.start < sreg.stop:
                if region.stop > sreg.stop:
                    sreg.stop = region.stop
                    modified = True
            elif region.stop > sreg.stop > region.start:
                sreg.stop = region.stop
                modified = True
        if not modified:
            current_regions.add(region)
        self.regions[region.chrom] = list(current_regions)

    def get_regions(self):
        region_list = []
        for chromlist in self.regions.values():
            region_list.extend(chromlist)
        return sorted(region_list)

    def get_inverse(self):
        inverse = defaultdict(tuple)
        chrom_names = get_grch38_chroms()
        for chrom in self.regions.keys():
            chrom_regions = self.regions[chrom]
            if len(chrom_regions) > 0:
                inv_reg = []
                next_start = 1
                for i, reg in enumerate(chrom_regions):
                    if i == 0 and int(reg.start) > 1:
                        inv_reg.append(VCFRegion(chrom, 1, reg.start))
                    elif int(reg.start) <= 1:
                        next_start = int(reg.stop)
                    else:
                        inv_reg.append(VCFRegion(chrom, next_start, reg.start))
                    next_start = int(reg.stop)
                if int(chrom_regions[-1].stop) < chrom_names[chrom]:
                    inv_reg.append(VCFRegion(chrom, next_start, chrom_names[chrom]))
                inverse[chrom] = tuple(inv_reg)
        inverse_list = []
        for chromlist in inverse.values():
            inverse_list.extend(chromlist)
        return sorted(inverse_list)
