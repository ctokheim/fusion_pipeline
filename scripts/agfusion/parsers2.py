import re
from agfusion.parsers import _Parser


class SimpleTSV(_Parser):
    def __init__(self, infile, logger):
        super(SimpleTSV, self).__init__(logger)

        fin = open(infile, 'r')
        for line in fin.readlines():

            line = line.strip().split('\t')
            if len(line) <= 2:
                continue

            self.fusions.append(
                {
                    'gene5prime': line[0],
                    'gene3prime': line[2],
                    'gene5prime_junction': int(line[1]),
                    'gene3prime_junction': int(line[3])
                }
            )
        fin.close()

        self._check_data()
