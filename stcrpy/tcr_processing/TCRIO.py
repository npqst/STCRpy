from Bio import PDB
from Bio.PDB.PDBIO import PDBIO


class TCRIO(PDBIO):
    def __init__(self):
        self.io = PDBIO()

    def save(
        self,
        tcr: "TCR",
        save_as: str = None,
        tcr_only: bool = False,
        format: str = "pdb",
    ):
        assert (
            tcr.__module__.split(".")[-1] == "TCR"
        ), f"{tcr} must be type TCR, not TCRStructure"
        structure_to_save = PDB.Model.Model(0)
        for chain in tcr.get_chains():
            chain.serial_num = 0
            structure_to_save.add(chain)
        if not tcr_only:
            for mhc in tcr.get_MHC():
                if (
                    mhc.__class__.__name__ == "MHCchain"
                ):  # handle MHC that is single chain
                    chain = mhc
                    chain.serial_num = 0
                    structure_to_save.add(chain)
                else:
                    for chain in mhc.get_chains():  # handle multichain MHC object
                        chain.serial_num = 0
                        structure_to_save.add(chain)

            for chain in tcr.get_antigen():
                chain.serial_num = 0
                structure_to_save.add(chain)

        self.io.set_structure(structure_to_save)
        if not save_as:
            if not tcr_only:
                save_as = f"{tcr.parent.parent.id}_{tcr.id}.{format}"
            else:
                save_as = f"{tcr.parent.parent.id}_{tcr.id}_TCR_only.{format}"

        self.io.save(save_as)
