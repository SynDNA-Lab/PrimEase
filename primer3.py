import re
import shutil
import pandas as pd 
from subprocess import Popen, PIPE
from dataclasses import dataclass, field

from selection import Selection
import subprocess

PRIMER_FASTA_FILE = "potential_primers.fasta"

@dataclass
class Primer3Interface:
    regions: list[str]
    primer3_path: str 
    result_path: str = field(default="primer3_result")
    primer_sites: pd.DataFrame = field(init=False)
    is_specific_file : bool = field(default=False)
    size_pcr_product : list[int] = field(default_factory=list)


    def run(self) -> None:
        self.write_settings()
        self.run_primer3()
        self.parse_results()
        self.generate_fasta()


    def write_settings(self) -> None:
        if self.is_specific_file == False : 
            # create a primer3 settings file
            shutil.copyfile("settings.bak", "settings")
            with open("settings", "a") as file:
                print(self.regions)
                file.write(str(self.regions)) # using __repr__ of ROIs
            with open("settings", "r") as f:
                lines = f.readlines()

            with open("settings", "w") as f:
                for line in lines:
                    if line.startswith("PRIMER_PRODUCT_SIZE_RANGE="):
                        f.write(f"PRIMER_PRODUCT_SIZE_RANGE={f"{self.size_pcr_product[0]}-{self.size_pcr_product[1]}"}\n")
                    else:
                        f.write(line)
            
        else : 
            shutil.copyfile("settings_spec.bak", "settings")
            with open("settings", "a") as file:
                file.write(str(self.regions)) # using __repr__ of ROIs
            with open("settings", "r") as f:
                lines = f.readlines()

            with open("settings", "w") as f:
                for line in lines:
                    if line.startswith("PRIMER_PRODUCT_SIZE_RANGE="):
                        f.write(f"PRIMER_PRODUCT_SIZE_RANGE={f"{self.size_pcr_product[0]}-{self.size_pcr_product[1]}"}\n")
                    else:
                        f.write(line)

    def run_primer3(self) -> None:
        cmd = f"{self.primer3_path} settings > {self.result_path}"
        process = Popen(args=cmd, stdout=PIPE, stderr=PIPE, shell=True)
        _, stderr = process.communicate()
        if stderr:
            print(stderr.decode("ascii"))


    def parse_results(self) -> None:
        # primer3 txt result to pandas dataframe
        cols = ["primer_pair_penalty",
                "primer_left_penalty",
                "primer_right_penalty",
                "primer_left_sequence",
                "primer_right_sequence",
                "primer_left",
                "primer_right",
                "primer_left_tm",
                "primer_right_tm",
                "primer_left_gc_percent",
                "primer_right_gc_percent",
                "primer_left_self_any_th",
                "primer_right_self_any_th",
                "primer_left_self_end_th",
                "primer_right_self_end_th",
                "primer_left_hairpin_th",
                "primer_right_hairpin_th",
                "primer_left_end_stability",
                "primer_right_end_stability",
                "primer_pair_compl_any_th",
                "primer_pair_compl_end_th",
                "primer_pair_product_size",
                "primer_pair_product_tm"]

        with open(self.result_path) as file:
            data = file.read()

        pattern = r"(?s)(?=PRIMER_PAIR_([0-9]+)_PENALTY)(.*?)(.*?PRIMER_PAIR_[0-9]+_PRODUCT_TM=[0-9]+.[0-9]+)"
        re_result = re.findall(pattern , data)

        if not re_result:
            raise ValueError(
                "Primer3 did not return any primer pairs. "
                "This may be due to overly strict design parameters, "
                "an unsuitable input sequence, or wrong primer3 configuration."
            )
        dfdat = []
        for r in re_result:
            d = r[-1].split("\n")
            e = [x for x in d if x]
            dfdat.append([x.split("=")[1] for x in e])

        df = pd.DataFrame(data=dfdat, columns=cols)
        
        names = []
        seqids = re.findall(r"SEQUENCE_ID=.+?(?=\n)", data)
        res =  re.findall(r"PRIMER_PAIR_NUM_RETURNED=.+?(?=\n)", data)

        for sid, rs in zip(seqids, res):
            s = sid.split("=")[1]
            for ids in range(int(rs.split("=")[1])):
                names.append(f"{s}_{ids}")

        positions = [nam.split("_")[-2] for nam in names]     
        df.insert(0, "name", names)
        df.insert(1, "position", positions)
        self.primer_sites = df


    def generate_fasta(self) -> str:
        buffer = ""
        for _, col in self.primer_sites.iterrows():
            buffer += f">{col['name']}_fwd\n{col['primer_left_sequence']}\n{col['name']}_rev\n{col['primer_right_sequence']}\n"

        with open(PRIMER_FASTA_FILE, "w") as file:
            file.write(buffer)
