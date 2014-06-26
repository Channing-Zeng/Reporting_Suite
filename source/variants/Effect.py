class Effect:
    """ EFF = Effect (
            Effect_Impact | Functional_Class | Codon_Change |
            Amino_Acid_Change | Amino_Acid_Length | Gene_Name |
            Transcript_BioType | Gene_Coding | Transcript_ID |
            Exon_Rank  | Genotype_Number [ | ERRORS | WARNINGS ]
        )
    """
    def __init__(self, string):
        self.efftype, rest = string[:-1].split('(')

        vals = rest.split('|')
        self.impact = vals[0]
        self.funclas = vals[1]
        self.cc = vals[2]
        self.aac = vals[3]
        self.pos = int(''.join(c for c in self.aac if c.isdigit())) if self.aac != '' else None
        self.aal = int(vals[4]) if vals[4] != '' else None
        self.gene = vals[5]
        self.biotype = vals[6]
        self.coding = vals[7] == 'CODING'
