from rdkit import Chem
from rdkit.Chem import AllChem
import itertools

DNA_BASES = u"ATGC"

def get_att_points(mol):
    att_points = []
    for a in mol.GetAtoms():
        if a.GetSymbol() == '*':
            att_points.append(a.GetIdx())

    return att_points
def map_idx(idx, idx_list, mol):
    abs_id = idx_list[idx]
    neigh_idx = mol.GetAtomWithIdx(abs_id).GetNeighbors()[0].GetIdx()
    return neigh_idx


class DNA_2_PDB(object):
    def __init__(self, ATCG_MOL):
        self.attach_point = Chem.MolFromSmiles('*')
        self.Na = Chem.MolFromSmiles('[Na+]')
        self.K = Chem.MolFromSmiles('[K+]')
        self.mol = ''

        self.atcg = ATCG_MOL
        self.atcg_att = [get_att_points(m) for _, m in enumerate(ATCG_MOL)]
        self.dna_list = ['A', 'T', 'C', 'G']

    def _add_motif(self, ac):

        cur_mol = Chem.ReplaceSubstructs(self.mol, self.attach_point, self.Na)[ac[0]]
        motif = self.atcg[ac[1]]
        att_point = self.atcg_att[ac[1]]
        # motif_atom = map_idx(ac[2], att_point, motif)
        motif = Chem.ReplaceSubstructs(motif, self.attach_point, self.K)[ac[2]]
        for atom in motif.GetAtoms():
            if atom.GetSymbol() == 'K':
                motif_atom = map_idx(0, [atom.GetIdx()], motif)
        motif = Chem.DeleteSubstructs(motif, self.K)
        next_mol = Chem.ReplaceSubstructs(cur_mol, self.Na, motif, replacementConnectionPoint=motif_atom)[0]
        self.mol = next_mol

    def mol_gen(self, DNA_Seqs):
        for DNA_Seq in DNA_Seqs:
            self.mol = self.atcg[self.dna_list.index(DNA_Seq[0])]
            for nuc in DNA_Seq[1:]:
                ac = [1, self.dna_list.index(nuc), 0]
                self._add_motif(ac)
            mol = Chem.DeleteSubstructs(self.mol, Chem.MolFromSmiles("*"))
            mol = AllChem.AddHs(mol)
            AllChem.EmbedMolecule(mol)
            # AllChem.MMFFOptimizeMolecule(mol)
            Chem.MolToMolFile(mol, '{}.mol'.format(DNA_Seq))

ATCG_SMI = ['*O[C@H]1C[C@H](N2C=NC3=C(N)N=CN=C32)O[C@@H]1COP(*)(=O)O',
            '*O[C@H]1C[C@H](N2C=C(C)C(=O)NC2=O)O[C@@H]1COP(*)(=O)O',
            '*O[C@H]1C[C@H](N2C=CC(N)=NC2=O)O[C@@H]1COP(*)(=O)O',
            '*O[C@H]1C[C@H](N2C=NC3=C2N=C(N)NC3=O)O[C@@H]1COP(*)(=O)O']
DNA_LIST = ['A', 'T', 'C', 'G']

ATCG_MOL = []
for smi in ATCG_SMI:
    ATCG_MOL.append(Chem.MolFromSmiles(smi))


dna2pdb = DNA_2_PDB(ATCG_MOL)
DNA_Seqs = ['ATGT']
# DNA_Seqs.extend(''.join(s) for s in itertools.product(DNA_BASES, repeat=5))
dna2pdb.mol_gen(DNA_Seqs)
