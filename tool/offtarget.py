import logging
import itertools
from dataclasses import dataclass, field

import pandas as pd 

from bowtie import BowtieResult
from config import Config



@dataclass
class OfftargetChecker:
    primer_candidates: pd.DataFrame
    bowtie_target: BowtieResult  
    bowtie_host: BowtieResult
    config: Config
    offtarget_ids: list[int] = field(init=False)
    final_primers : pd.DataFrame = field(init=False)


    def __post_init__(self) -> None:
        delete_buffer = []
        delete_buffer.extend(self.get_sponges(self.bowtie_target.result))

        delete_buffer.extend(self.filter_offtarget(self.bowtie_target.result))
        
        self.offtarget_ids = list(set(delete_buffer))
        logging.info(f"{len(self.offtarget_ids)} Offtargets found")

        host_ids = self.get_host_hits(self.bowtie_host.result)
        logging.info(f"{len(host_ids)} primers align to host genome (will be removed)")
        delete_buffer.extend(host_ids)

        self.create_primer_list()


    def get_sponges(self, data_frame:pd.DataFrame) -> list[str]:
        grouped = data_frame.groupby(["id", "orientation"])
        sponge_df = grouped.filter(lambda x: len(x)>= self.config.sponge_value)
        return sponge_df.id.unique()


    def filter_offtarget(self, data_frame:pd.DataFrame) -> list[str]:
        grouped = data_frame.groupby(["id", "reference"])
        df_filtered = grouped.filter(lambda x: len(x)>=2)
        return_list = []
        for _, group in df_filtered.groupby(["id", "reference"]):
            if len(group.strand.unique())<2:
                continue
            fwd = group[group.strand=="-"].start.unique()
            rev = group[group.strand=="+"].start.unique()
            
            for f, r in itertools.product(fwd, rev):
                difference = r-f
                if (difference > 0) and (difference <= self.config.offtarget_size_cutoff):
                    return_list.append(group.id.values[0])

        return return_list

    def get_host_hits(self, data_frame: pd.DataFrame) -> list[str]:
        """Return all primer IDs that align anywhere to the host genome."""
        if "id" not in data_frame.columns:
            raise ValueError("bowtie_host.result must contain an 'id' column.")
        return data_frame.id.unique().tolist()

    def create_primer_list(self) -> None:
        offtargets = [f"target_{ids}" for ids in self.offtarget_ids]
        final_df = self.primer_candidates[~self.primer_candidates.name.isin(offtargets)]
        final_df = final_df.groupby(["position"]).head(self.config.top)
        final_df.to_csv("qTagGer_Output.csv", index=False)
        missed_primers = set(self.primer_candidates.position.unique()) - set(final_df.position.unique())
        [print(f"Unable to find primer for position {mp}") for mp in missed_primers]
        self.final_primers = final_df
            
