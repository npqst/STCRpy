# Adapted from https://github.com/bjornwallner/DockQ

# MIT License applies:
# MIT License

# Copyright (c) 2024 bjornwallner

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import warnings
from collections import namedtuple
import itertools
import Bio
import json
from functools import partial
import logging
import numpy as np
import copy
import os


class TCRDockQ:
    def __init__(self, TCR_pMHC_interface=True, **kwargs):
        """
        Initialize the TCRDockQ class. See `Dockq --help` for more information.

        Args:
            TCR_pMHC_interface (bool): Whether to use evaluate the entire TCR-pMHC interface or evaluate spearate chains. Deaults to True. This will create a merged psudo-chain for the antigen and the TCR.
            **kwargs: Keyword arguments from DockQ.
        """

        DockQArgs = namedtuple(
            "DockQArgs",
            [
                "mapping",
                "small_molecule",
                "capri_peptide",
                "no_align",
                "json",
                "n_cpu",
                "max_chunk",
                "allowed_mismatches",
                "verbose",
                "short",
                "print_results",
            ],
        )
        self.TCR_pMHC_interface = TCR_pMHC_interface
        self.args = DockQArgs(
            mapping=kwargs.get("mapping", None),
            small_molecule=kwargs.get("small_molecule", None),
            capri_peptide=kwargs.get("capri_peptide", False),
            no_align=kwargs.get("no_align", False),
            json=kwargs.get("json", None),
            n_cpu=kwargs.get("n_cpu", 8),
            max_chunk=kwargs.get("max_chunk", 512),
            allowed_mismatches=kwargs.get("allowed_mismatches", 0),
            verbose=kwargs.get("verbose", False),
            short=kwargs.get("short", False),
            print_results=kwargs.get("print_results", False),
        )
        return

    def try_tcr_mapping(self, dock: "abTCR", reference: "abTCR"):
        try:
            dock_chain_mapping = {v: k for k, v in dock.get_chain_mapping().items()}
            ref_chain_mapping = {v: k for k, v in reference.get_chain_mapping().items()}

            mapping = [
                (dock_chain_mapping[chain_type], ref_chain_mapping[chain_type])
                for chain_type in dock_chain_mapping
                if chain_type in ref_chain_mapping
            ]
            mapping = list(zip(*mapping))
            string_mapping = f"{''.join(mapping[0])}:{''.join(mapping[1])}"
            return mapping, string_mapping
        except Exception:
            return None, None

    def retrieve_chain(self, tcr: "TCR", chain_id: str):
        chain_map = tcr.get_chain_mapping()
        if chain_map[chain_id] in ["VA", "VB", "VD", "VG"]:
            return tcr[chain_id]
        elif chain_map[chain_id] == "Ag":
            return [ag for ag in tcr.get_antigen() if ag.id == chain_id][0]
        else:
            return tcr.get_MHC()[0][chain_id]

    def filter_residues(self, dock: "abTCR", reference: "abTCR", mapping):
        from stcrpy.tcr_formats.tcr_formats import get_sequences

        filtered_dock = copy.deepcopy(dock)
        filtered_reference = copy.deepcopy(reference)

        for i, c_id in enumerate(mapping[0]):
            dock_chain = self.retrieve_chain(dock, c_id)
            ref_chain = self.retrieve_chain(reference, mapping[1][i])

            # get the sequence of the dock chain and ref chain, then align them. Identify any mismathces in the sequence and remove those residues from the chain
            dock_seq = get_sequences(dock_chain, amino_acids_only=False)[dock_chain.id]
            ref_seq = get_sequences(ref_chain, amino_acids_only=False)[ref_chain.id]
            aligner = Bio.Align.PairwiseAligner(match=1, mismatch=-1, gap_score=-1)
            alignment = aligner.align(dock_seq, ref_seq)[0]
            dock_indices, ref_indices = alignment.indices
            filtered_dock_residues = [
                dock_chain.child_list[idx].copy()
                for idx in dock_indices[
                    np.logical_and(ref_indices != -1, dock_indices != -1)
                ]
            ]
            filtered_ref_residues = [
                ref_chain.child_list[idx].copy()
                for idx in ref_indices[
                    np.logical_and(ref_indices != -1, dock_indices != -1)
                ]
            ]

            filtered_dock_chain = self.retrieve_chain(filtered_dock, c_id)
            filtered_ref_chain = self.retrieve_chain(filtered_reference, mapping[1][i])

            # remove residues
            for r in dock_chain:
                filtered_dock_chain.detach_child(r.id)
            for r in ref_chain:
                filtered_ref_chain.detach_child(r.id)

            # add filtered residues that aligned in dock and ref
            for r in filtered_dock_residues:
                filtered_dock_chain.add(r)
            for r in filtered_ref_residues:
                filtered_ref_chain.add(r)

        return filtered_dock, filtered_reference

    def tcr_to_structure(
        self, tcr: "abTCR", mapping_filter=None, merge_chains=True
    ) -> Bio.PDB.Structure.Structure:
        # structure = Bio.PDB.Structure.Structure(tcr.parent.parent.id)
        from stcrpy import tcr_formats

        model = Bio.PDB.Model.Model(tcr.parent.id)

        if merge_chains:
            tcr_chains_to_merge = [
                c for c in tcr.get_chains() if c.id in mapping_filter
            ]
            tcr_chain = tcr_formats.tcr_formats.merge_chains(
                tcr_chains_to_merge, new_chain_id="T"
            )
            tcr_chain.is_het = False
            tcr_chain.sequence = tcr_formats.tcr_formats.get_sequences(tcr_chain)["T"]
            model.add(tcr_chain)

            antigen_MHC_chains = tcr.get_antigen()
            antigen_MHC_chains.extend(
                [c for m in tcr.get_MHC() for c in m.get_chains()]
            )
            antigen_MHC_chains_to_merge = [
                c for c in antigen_MHC_chains if c.id in mapping_filter
            ]
            pMHC_chain = tcr_formats.tcr_formats.merge_chains(
                antigen_MHC_chains_to_merge, new_chain_id="M"
            )
            pMHC_chain.is_het = False
            pMHC_chain.sequence = tcr_formats.tcr_formats.get_sequences(pMHC_chain)["M"]
            model.add(pMHC_chain)

        else:
            sequences = tcr_formats.tcr_formats.get_sequences(tcr)
            for chain in tcr.get_chains():
                chain.is_het = False
                chain.sequence = sequences[chain.id]
                model.add(chain)
            for antigen in tcr.get_antigen():
                sequences = tcr_formats.tcr_formats.get_sequences(antigen)
                antigen.is_het = False
                antigen.sequence = sequences[antigen.id]
                model.add(antigen)
            for mhc in tcr.get_MHC():
                sequences = tcr_formats.tcr_formats.get_sequences(mhc)
                for chain in mhc.get_chains():
                    chain.is_het = False
                    chain.sequence = sequences[chain.id]
                    model.add(chain)
        # structure.add(model)
        return model

    def tcr_dockq(self, dock: "abTCR", reference: "abTCR", save_merged_complex: bool=False) -> float:
        """
        Calculate DockQ metrics for a TCR-pMHC complex.

        This method evaluates the quality of a predicted TCR-pMHC complex (the "dock" structure)
        against a reference (native) structure using the DockQ scoring system. It supports both
        merged TCR-pMHC interfaces and separate chain evaluations, depending on the class setting.

        Args:
            dock (abTCR): The predicted TCR-pMHC complex structure to be evaluated.
            reference (abTCR): The reference (native) TCR-pMHC complex structure.
            save_merged_complex (bool, optional): If True, saves the merged complex structures
                as PDB files for further inspection. Defaults to False.

        Returns:
            dict: A dictionary containing DockQ results, including DockQ score, F1, iRMSD, LRMSD,
                fnat, fnonnat, and other relevant metrics for the best mapping.
        """

        try:
            import DockQ
            from DockQ import DockQ as dockq
        except (ImportError, ModuleNotFoundError) as e:
            warnings.warn(
                f"DockQ is not installed. Error: {e}. Please install it using `pip install DockQ`."
            )

        mapping, string_mapping = self.try_tcr_mapping(dock, reference)

        filtered_dock, filtered_reference = self.filter_residues(
            dock, reference, mapping
        )

        model_structure = self.tcr_to_structure(
            filtered_dock,
            mapping_filter=mapping[0],
            merge_chains=self.TCR_pMHC_interface,
        )
        native_structure = self.tcr_to_structure(
            filtered_reference,
            mapping_filter=mapping[1],
            merge_chains=self.TCR_pMHC_interface,
        )

        # save the structures to pdb files
        import Bio

        io = Bio.PDB.PDBIO()
        io.set_structure(model_structure)
        io.save(f"model_structure_{dock.parent.parent.id}_{dock.id}.pdb")
        io.set_structure(native_structure)
        io.save(f"native_structure_{reference.parent.parent.id}_{reference.id}.pdb")

        if self.TCR_pMHC_interface:
            string_mapping = "TM:TM"

        initial_mapping, model_chains, native_chains = dockq.format_mapping(
            string_mapping, small_molecule=self.args.small_molecule
        )

        model_structure = dockq.load_PDB(
            f"model_structure_{dock.parent.parent.id}_{dock.id}.pdb",
            chains=model_chains,
            small_molecule=self.args.small_molecule,
        )
        native_structure = dockq.load_PDB(
            f"native_structure_{reference.parent.parent.id}_{reference.id}.pdb",
            chains=native_chains,
            small_molecule=self.args.small_molecule,
        )

        if not save_merged_complex:
            os.remove(f"model_structure_{dock.parent.parent.id}_{dock.id}.pdb")
            os.remove(f"native_structure_{reference.parent.parent.id}_{reference.id}.pdb")

        # check user-given chains are in the structures
        model_chains = (
            [c.id for c in model_structure] if not model_chains else model_chains
        )
        native_chains = (
            [c.id for c in native_structure] if not native_chains else native_chains
        )

        if len(model_chains) < 2 or len(native_chains) < 2:
            print("Need at least two chains in the two inputs\n")
            return

        # permute chains and run on a for loop
        best_dockq = -1
        best_result = None
        best_mapping = None

        model_chains_to_combo = [
            mc for mc in model_chains if mc not in initial_mapping.values()
        ]
        native_chains_to_combo = [
            nc for nc in native_chains if nc not in initial_mapping.keys()
        ]

        chain_clusters, reverse_map = dockq.group_chains(
            model_structure,
            native_structure,
            model_chains_to_combo,
            native_chains_to_combo,
            self.args.allowed_mismatches,
        )
        chain_maps = dockq.get_all_chain_maps(
            chain_clusters,
            initial_mapping,
            reverse_map,
            model_chains_to_combo,
            native_chains_to_combo,
        )

        num_chain_combinations = dockq.count_chain_combinations(chain_clusters)
        # copy iterator to use later
        chain_maps, chain_maps_ = itertools.tee(chain_maps)

        low_memory = num_chain_combinations > 100
        run_chain_map = partial(
            dockq.run_on_all_native_interfaces,
            model_structure,
            native_structure,
            no_align=self.args.no_align,
            capri_peptide=self.args.capri_peptide,
            low_memory=low_memory,
        )

        if num_chain_combinations > 1:
            cpus = min(num_chain_combinations, self.args.n_cpu)
            chunk_size = min(
                self.args.max_chunk, max(1, num_chain_combinations // cpus)
            )

            # for large num_chain_combinations it should be possible to divide the chain_maps in chunks
            result_this_mappings = dockq.progress_map(
                run_chain_map,
                chain_maps,
                total=num_chain_combinations,
                n_cpu=cpus,
                chunk_size=chunk_size,
            )

            for chain_map, (result_this_mapping, total_dockq) in zip(
                chain_maps_, result_this_mappings
            ):

                if total_dockq > best_dockq:
                    best_dockq = total_dockq
                    best_result = result_this_mapping
                    best_mapping = chain_map

            if (
                low_memory
            ):  # retrieve the full output by rerunning the best chain mapping
                best_result, total_dockq = dockq.run_on_all_native_interfaces(
                    model_structure,
                    native_structure,
                    chain_map=best_mapping,
                    no_align=self.args.no_align,
                    capri_peptide=self.args.capri_peptide,
                    low_memory=False,
                )

        else:  # skip multi-threading for single jobs (skip the bar basically)
            best_mapping = next(chain_maps)
            best_result, best_dockq = run_chain_map(best_mapping)

        if not best_result:
            logging.error(
                "Could not find interfaces in the native model. Please double check the inputs or select different chains with the --mapping flag."
            )
            return

        info = dict()
        info["model"] = f"{dock.parent.parent.id}_{dock.id}"
        info["native"] = f"{reference.parent.parent.id}_{reference.id}"
        info["best_dockq"] = best_dockq
        info["best_result"] = best_result
        info["GlobalDockQ"] = best_dockq / len(best_result)
        info["best_mapping"] = best_mapping
        info["best_mapping_str"] = f"{dockq.format_mapping_string(best_mapping)}"

        if self.args.json:
            if not isinstance(self.args.json, str):
                json_file = f"dockq_{model_structure.id}_{dock.id}_{native_structure.id}_{reference.id}.json"
            else:
                json_file = self.args.json
            with open(json_file, "w") as fp:
                json.dump(info, fp)

        if self.args.print_results:
            dockq.print_results(
                info,
                self.args.short,
                self.args.verbose,
                self.args.capri_peptide,
                self.args.small_molecule,
            )

        return info
