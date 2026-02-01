import pandas as pd 
from subprocess import Popen, PIPE
from dataclasses import dataclass, field

from config import Config
import os


BASENAME = "target"

@dataclass
class BowtieResult:
    result_path: str
    result: pd.DataFrame = field(init=False)

    def __post_init__(self) -> None:
        self.parse_data()


    def split_data(self, lst:list[str], instruction:str) -> list[str]:
        if instruction == "orientation":
            return [l.split("_")[-1] for l in lst]
        elif instruction == "id":
            return ["_".join([l.split("_")[-3], l.split("_")[-2]]) for l in lst]
        else:
            raise ValueError(f"Unknown bowtie parser instruction {instruction}")

    def split_data_corrected(self, lst: list[str], instruction: str) -> list[str]:
        if instruction == "orientation":
            return [l.split("_")[-1] if len(l.split("_")) >= 2 else None for l in lst]
        elif instruction == "id":
            return [
                "_".join([l.split("_")[-3], l.split("_")[-2]])
                if len(l.split("_")) >= 3 else None
                for l in lst
            ]
        else:
            raise ValueError(f"Unknown bowtie parser instruction {instruction}")


    def parse_data(self) -> None:
        expected_columns = ["name","strand","reference","start","sequence","quality","instances"]
        
        if not os.path.exists(self.result_path) or os.stat(self.result_path).st_size == 0:
            # empty dataframe but include id/orientation
            df = pd.DataFrame(columns=expected_columns + ["id","orientation"])
        else:
            # read all columns safely
            df = pd.read_csv(self.result_path, sep="\t", header=None, dtype=str)
            # keep only first 7 columns
            df = df.iloc[:, :7] if df.shape[1] >= 7 else pd.DataFrame(columns=expected_columns)
            df.columns = expected_columns
            df = df.drop_duplicates()
        
            # always add id/orientation even if df is empty
            df["id"] = df["name"].apply(lambda x: "_".join(x.split("_")[-3:-1]) if pd.notnull(x) else None)
            df["orientation"] = df["name"].apply(lambda x: x.split("_")[-1] if pd.notnull(x) else None)

        self.result = df







@dataclass
class BowtieInterface:
    config: Config
    random_identifiers : list[str] = field(default_factory=list)
    output_target: str = field(default="bt_target.csv")
    output_host : str = field(default ="bt_host.csv")
    result_target: BowtieResult = field(init=False)
    result_host: BowtieResult  = field(init=False)


    def __post_init__(self) -> None:
        self.create_index()

        self.run_bowtie(
            index=f"{self.config.bowtie_path}/target/{BASENAME}", 
            fasta_path="potential_primers.fasta", 
            output_path=self.output_target)
        host_dfs = []
        for genome_id in self.random_identifiers:
            host_index = f"{self.config.bowtie_path}/host/{genome_id}"
            host_output_path = f"bt_host_{genome_id}.csv"
            if os.path.exists(f"bowtie_index/host/{genome_id}.1.ebwt"):
                self.run_bowtie(index=host_index, fasta_path="potential_primers.fasta", output_path=host_output_path)
                host_dfs.append(pd.read_csv(host_output_path, sep="\t", header=None))

        # Concatenate
        print(host_dfs)
        if host_dfs:
            combined_df = pd.concat(host_dfs, ignore_index=True)
            combined_df.to_csv(self.output_host, index=False, header=False)
        else:
            # empty DataFrame if no host genomes
            combined_df = pd.DataFrame()
            combined_df.to_csv(self.output_host, index=False)
        # Store
        self.result_host = BowtieResult(result_path=self.output_host)
        self.result_target = BowtieResult(result_path=self.output_target)
        
    

    def run_command(self, cmd: str) -> None:
        process = Popen(args=cmd, stdout=PIPE, stderr=PIPE, shell=True)
        _, stderr = process.communicate()
        if stderr:
            print(stderr.decode("ascii"))
    

    def create_index(self) -> None:
        btpath = self.config.bowtie_path 
        targetpath = "../../target.fasta" 

        cmd = f"cd {btpath}/target && bowtie-build -f {targetpath} {BASENAME} && cd ../.."
        self.run_command(cmd=cmd)
        
        for genome_id in self.random_identifiers:
            host_fasta = f"bowtie_index/host/{genome_id}.fasta"
            host_index_file = f"bowtie_index/host/{genome_id}.1.ebwt"
            if os.path.exists(host_fasta) and not os.path.exists(host_index_file):
                cmd = f"cd {btpath}/host && bowtie-build -f {genome_id}.fasta {genome_id} && cd ../.."
                self.run_command(cmd=cmd)
        
    def run_bowtie(self, index:str, fasta_path:str, output_path:str) -> None:
        cmd = f"bowtie -x {index} -a -f {fasta_path} -v 3 > {output_path}"
        self.run_command(cmd=cmd)
