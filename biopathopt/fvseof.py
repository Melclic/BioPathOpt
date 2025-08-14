#Taken and modified from https://github.com/LucasCoppens/fvseof

from __future__ import annotations


import logging
from typing import Dict, List

import cobra
from cobra import Reaction, Metabolite
import numpy as np
import pandas as pd
from cobra import Reaction
from cobra.flux_analysis import flux_variability_analysis
from cobra.flux_analysis.loopless import loopless_solution
from tqdm import tqdm


class FVSEOF:
    """Flux Variability Scanning based on Enforced Objective Flux (FVSEOF).

    This class implements the FVSEOF algorithm to identify reactions whose
    fluxes change as the production of a target metabolite is gradually enforced.

    Attributes:
        model (cobra.Model): A copy of the input model used internally.
        biomass_reaction_id (str): Reaction ID of the biomass objective.
        target_metabolite_id (str): Metabolite ID of the desired product.
        product_sink_reaction_id (str): ID of the added sink reaction for product.
        product_max_theoretical_yield (float): Maximal theoretical product flux.
        maximal_biomass_growth (float): Unconstrained maximal growth rate.
        check_essential_reaction_threshold (float): Threshold fraction of
            maximal growth below which a reaction knockout is considered essential.
    """

    def __init__(
        self,
        model: cobra.Model,
        biomass_reaction_id: str,
        target_metabolite_id: str,
        essential_reaction_threshold: float = 0.5,
    ) -> None:
        """Initialize the FVSEOF class.

        Args:
            model: COBRA model to analyze.
            biomass_reaction_id: Reaction ID for the biomass objective.
            target_metabolite_id: Metabolite ID for which a sink will be enforced.
            essential_reaction_threshold: Fraction of maximal biomass growth
                below which a knockout is considered essential. Defaults to 0.5.

        Raises:
            AssertionError: If the biomass reaction or target metabolite is missing.

        Example:
            >>> fvseof = FVSEOF(model, "BIOMASS_Ecoli_core_w_GAM", "etoh_c")
        """
        self.model: cobra.Model = model.copy()
        self.biomass_reaction_id: str = biomass_reaction_id
        self.target_metabolite_id: str = target_metabolite_id

        assert self.biomass_reaction_id in [r.id for r in self.model.reactions], (
            "Biomass reaction not in model."
        )
        assert self.target_metabolite_id in [m.id for m in self.model.metabolites], (
            "Target metabolite not in model."
        )

        self.product_sink_reaction_id: str = self.add_product_sink_reaction(
            self.target_metabolite_id
        )
        self.product_max_theoretical_yield: float = (
            self.calc_product_maximal_theoretical_yield()
        )
        self.maximal_biomass_growth: float = self.calc_maximal_biomass_growth()

        self.check_essential_reaction_threshold: float = essential_reaction_threshold



    def find_sink_reaction_ids(self, metabolite_id: str) -> List[str]:
        """Return IDs of all sink reactions that consume a given metabolite.

        A **sink reaction** here is defined structurally as a reaction that:
          1) Involves **exactly one** metabolite (the queried metabolite), and
          2) Is **irreversible** in the **consuming direction** (i.e., consumes the
             metabolite into nothing, without producing any other metabolites).

        We accept either orientation depending on stoichiometry and bounds:
          - If the stoichiometric coefficient for the metabolite is **negative**,
            then positive flux consumes it, so we require bounds **lb >= 0 and ub > 0**.
          - If the stoichiometric coefficient is **positive**, then negative flux
            consumes it, so we require bounds **lb < 0 and ub <= 0**.

        Args:
            model: COBRApy model to search.
            metabolite_id: Metabolite ID whose sink reactions are sought.

        Returns:
            List of reaction IDs that are sink reactions for the metabolite.

        Raises:
            KeyError: If `metabolite_id` is not present in the model.

        Example:
            >>> ids = find_sink_reaction_ids(model, "etoh_c")
            >>> ids
            ['etoh_c_sink']
        """
        # Ensure metabolite exists (raises KeyError if not found).
        _met: Metabolite = self.model.metabolites.get_by_id(metabolite_id)

        sink_ids: List[str] = []
        for rxn in self.model.reactions:  # type: Reaction
            mets = rxn.metabolites
            if len(mets) != 1:
                continue  # must involve exactly one metabolite

            only_met, coeff = next(iter(mets.items()))
            if only_met.id != metabolite_id:
                continue  # must be the queried metabolite

            # Skip blocked reactions
            lb, ub = rxn.lower_bound, rxn.upper_bound
            if lb == 0 and ub == 0:
                continue

            # Determine if bounds permit only consumption (irreversible in consuming direction)
            coeff = float(coeff)
            consumes_with_positive_flux = coeff < 0 and lb >= 0 and ub > 0
            consumes_with_negative_flux = coeff > 0 and lb < 0 and ub <= 0

            if consumes_with_positive_flux or consumes_with_negative_flux:
                sink_ids.append(rxn.id)

        return sink_ids


    def add_product_sink_reaction(self, target_metabolite_id: str) -> str:
        """Add a sink reaction for the target metabolite.

        Args:
            target_metabolite_id: Metabolite ID to drain via a sink.

        Returns:
            The ID of the created sink reaction.

        Example:
            >>> sink_id = fvseof.add_product_sink_reaction("etoh_c")
        """
        #check that there isn't already a sink for a given metabolite id
        reacs = self.find_sink_reaction_ids(target_metabolite_id)
        if len(reacs)==1:
            return reacs[0]
        elif len(reacs)>1:
            logging.warning(f'There are multiple sinks for {target_metabolite_id}: {reacs}')
        reaction_id = f"{target_metabolite_id}_fvseof_sink"
        product_sink_reaction: Reaction = Reaction(reaction_id)
        product_sink_reaction.add_metabolites(
            {self.model.metabolites.get_by_id(target_metabolite_id): -1}
        )
        product_sink_reaction.bounds = (0.0, 1000.0)
        self.model.add_reactions([product_sink_reaction])
        return reaction_id

    def calc_product_maximal_theoretical_yield(self) -> float:
        """Compute maximal theoretical production rate for the target product.

        Returns:
            Maximal theoretical product flux (1/h if model is in standard units).

        Example:
            >>> max_prod = fvseof.calc_product_maximal_theoretical_yield()
        """
        with self.model as model_ctx:
            model_ctx.objective = self.product_sink_reaction_id
            sol = loopless_solution(model_ctx)
            return float(sol.fluxes[self.product_sink_reaction_id])

    def calc_maximal_biomass_growth(self) -> float:
        """Compute maximal growth with unconstrained product formation.

        Returns:
            Maximal growth rate (1/h).

        Example:
            >>> mu_max = fvseof.calc_maximal_biomass_growth()
        """
        with self.model as model_ctx:
            model_ctx.objective = self.biomass_reaction_id
            sol = loopless_solution(model_ctx)
            return float(sol.fluxes[self.biomass_reaction_id])

    def check_essential_reaction(self, reaction_id: str) -> bool:
        """Determine if a reaction is essential for growth.

        The reaction is considered essential if knocking it out reduces the
        growth rate below `check_essential_reaction_threshold * maximal_biomass_growth`.

        Args:
            reaction_id: Reaction ID to test.

        Returns:
            True if essential, False otherwise. If optimization fails, returns True.

        Example:
            >>> is_essential = fvseof.check_essential_reaction("PFK")
        """
        with self.model as model_ctx:
            model_ctx.objective = self.biomass_reaction_id
            reaction = model_ctx.reactions.get_by_id(reaction_id)
            # Knockout by fixing bounds to zero
            reaction.lower_bound = 0.0
            reaction.upper_bound = 0.0

            try:
                sol = loopless_solution(model_ctx)
                mu = float(sol.fluxes[self.biomass_reaction_id])
                return (
                    mu / self.maximal_biomass_growth
                    < self.check_essential_reaction_threshold
                )
            except Exception:
                # If solver fails, be conservative and mark as essential
                return True

    def run(
        self,
        n_steps: int = 10,
        check_essentiality: bool = False,
        fva: bool = True,
        fva_n_processes: int = 1,
        use_progressbar: bool = False,
        fraction_of_optimum: float = 0.95,
    ) -> pd.DataFrame:
        """Run FVSEOF across increasing enforced product flux levels.

        For each step, the lower bound of the product sink reaction is increased
        linearly from 0 to the maximal theoretical yield (exclusive in the last step
        if using `range(0, n_steps)`). Fluxes are recorded via FVA (mean of min/max)
        or FBA (loopless solution). Reactions are classified by how their flux
        changes: Up, Down, or Reverse.

        Args:
            n_steps: Number of product enforcement steps. Defaults to 10.
            check_essentiality: If True, evaluate essentiality for targets. Defaults to False.
            fva: If True, use FVA (FVSEOF). If False, use FBA (FSEOF). Defaults to True.
            fva_n_processes: Parallel processes for FVA. Defaults to 1.
            use_progressbar: If True, show tqdm progress bars. Defaults to False.

        Returns:
            DataFrame indexed by reaction_id with columns:
                - target_type: {"Up","Down","Reverse"}
                - reaction: Human-readable reaction name
                - reaction formula: Reaction string with metabolite names
                - gene_reaction_rule: GPR string (OR->"/", AND->"+")
                - slope: Linear slope of flux vs. enforced product LB
                - step_i: Flux value at step i (mean for FVA or FBA flux)
                - essentiality (optional): True/False if `check_essentiality` is True

        Example:
            >>> df_targets = fvseof.run(n_steps=8, check_essentiality=True, fva=True)
            >>> df_targets.head()
        """
        steps_lower_bounds: List[float] = [
            (n / n_steps) * self.product_max_theoretical_yield for n in range(0, n_steps)
        ]
        method = "FVA" if fva else "FBA"

        # Collect per-step fluxes for all reactions
        per_step_fluxes: Dict[str, List[float]] = {
            r.id: [] for r in self.model.reactions
        }

        steps_iterator = (
            tqdm(enumerate(steps_lower_bounds), desc="FVSEOF", total=len(steps_lower_bounds))
            if use_progressbar
            else enumerate(steps_lower_bounds)
        )

        for i, step_lower_bound in steps_iterator:
            logging.info(
                "Performing %s for step %d/%d...",
                method,
                i + 1,
                len(steps_lower_bounds),
            )
            with self.model as model_ctx:
                model_ctx.objective = self.biomass_reaction_id
                model_ctx.reactions.get_by_id(
                    self.product_sink_reaction_id
                ).lower_bound = step_lower_bound

                if fva:
                    fva_df = flux_variability_analysis(
                        model_ctx, fraction_of_optimum=fraction_of_optimum, processes=fva_n_processes
                    )
                    for r in self.model.reactions:
                        per_step_fluxes[r.id].append(
                            (float(fva_df.loc[r.id, "maximum"]) + float(fva_df.loc[r.id, "minimum"])) / 2.0
                        )
                else:
                    fba_sol = loopless_solution(model_ctx)
                    for r in self.model.reactions:
                        per_step_fluxes[r.id].append(float(fba_sol.fluxes[r.id]))

        logging.info("FVSEOF sweep complete.")

        # Classify targets based on first vs last flux values
        target_types: Dict[str, str] = {}
        for r in self.model.reactions:
            r_id = r.id
            first = per_step_fluxes[r_id][0]
            last = per_step_fluxes[r_id][-1]
            if first == last:
                continue
            if first * last >= 0:
                target_types[r_id] = "Up" if abs(last) > abs(first) else "Down"
            else:
                target_types[r_id] = "Reverse"

        # Linear slope of flux vs enforced product flux
        slopes: Dict[str, float] = {}
        for r in self.model.reactions:
            r_id = r.id
            #TODO: when the fit is too poor then ignore
            slopes[r_id] = float(np.polyfit(steps_lower_bounds, per_step_fluxes[r_id], 1)[0])

        # Optional essentiality check
        essentialities: Dict[str, bool] = {}
        if check_essentiality:
            logging.info("Calculating essentialities...")
            reac_iter = (
                tqdm(self.model.reactions, desc="Essentiality")
                if use_progressbar
                else self.model.reactions
            )
            for r in reac_iter:
                if r.id in target_types:
                    essentialities[r.id] = self.check_essential_reaction(r.id)

        # Build result DataFrame
        df_data: Dict[str, Dict[str, float | str | bool]] = {}
        for r in self.model.reactions:
            r_id = r.id
            if r_id not in target_types:
                continue

            entry: Dict[str, float | str | bool] = {
                "target_type": target_types[r_id],
                "reaction": r.name,
                "reaction formula": r.build_reaction_string(use_metabolite_names=True),
                "gene_reaction_rule": r.gene_reaction_rule.replace(" or ", "/").replace(" and ", "+"),
                "slope": slopes[r_id],
            }

            if check_essentiality:
                entry["essentiality"] = essentialities.get(r_id, False)

            for i, _lb in enumerate(steps_lower_bounds):
                entry[f"step_{i}"] = per_step_fluxes[r_id][i]

            df_data[r_id] = entry

        df = pd.DataFrame.from_dict(df_data, orient="index")
        df.index.name = "reaction_id"

        # Sort: target_type desc (lexicographically puts Up first), slope by abs desc
        df = df.sort_values(
            by=["target_type", "slope"],
            key=lambda col: col.abs() if col.name == "slope" else col,
            ascending=[False, False],
        )
        return df
