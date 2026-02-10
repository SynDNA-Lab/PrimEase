import pandas as pd 
from subprocess import Popen, PIPE
from dataclasses import dataclass, field

from config import Config
import shutil
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
        expected_columns = [
            "name", "strand", "reference",
            "start", "sequence", "quality", "instances"
        ]

        # File missing or truly empty
        if not os.path.exists(self.result_path) or os.stat(self.result_path).st_size == 0:
            self.result = pd.DataFrame(columns=expected_columns + ["id", "orientation"])
            return

        try:
            df = pd.read_csv(self.result_path, sep="\t", header=None, dtype=str)
        except pd.errors.EmptyDataError:
            self.result = pd.DataFrame(columns=expected_columns + ["id", "orientation"])
            return

        # Not enough columns
        if df.shape[1] < 7:
            self.result = pd.DataFrame(columns=expected_columns + ["id", "orientation"])
            return

        df = df.iloc[:, :7]
        df.columns = expected_columns
        df = df.drop_duplicates()

        df["start"] = pd.to_numeric(df["start"], errors="coerce")

        df["id"] = df["name"].apply(
            lambda x: "_".join(x.split("_")[-3:-1]) if pd.notnull(x) else None
        )
        df["orientation"] = df["name"].apply(
            lambda x: x.split("_")[-1] if pd.notnull(x) else None
        )

        self.result = df








@dataclass
class BowtieInterface:
    config: Config
    task_dir : str = field(default="None")
    target_dir : str = field(init = False)
    host_dir : str = field(init = False)
    output_target: str = field(default="bt_target.csv")
    output_host : str = field(default ="bt_host.csv")
    result_target: BowtieResult = field(init=False)
    result_host: BowtieResult  = field(init=False)


    def __post_init__(self) -> None:
        Index = self.create_index()

        self.run_bowtie(
            index=f"{self.target_dir}/target", 
            fasta_path="potential_primers.fasta", 
            output_path=self.output_target)
        
        host_dfs = []

        for index_elem in Index:

            host_index = os.path.join(self.task_dir, f"host/{index_elem}")
            host_output_path = f"bt_host_{index_elem}.csv"
            print(host_index)
            if os.path.exists(f"{host_index}.1.ebwt"):
                self.run_bowtie(index=host_index, fasta_path="potential_primers.fasta", output_path=host_output_path)
                if os.path.exists(host_output_path) and os.stat(host_output_path).st_size > 0:
                    host_dfs.append(
                        pd.read_csv(host_output_path, sep="\t", header=None)
                        )
                if os.path.exists(host_output_path):
                    os.remove(host_output_path)
                    

        # Concatenate
        if host_dfs:
            combined_df = pd.concat(host_dfs, ignore_index=True)
            combined_df.to_csv(self.output_host, index=False, header=False, sep="\t")
        else:
            # empty DataFrame if no host genomes
            combined_df = pd.DataFrame()
            combined_df.to_csv(self.output_host, index=False)
        self.result_host = BowtieResult(result_path=self.output_host)
        self.result_target = BowtieResult(result_path=self.output_target)
        
    

    def run_command(self, cmd: str) -> None:
        process = Popen(args=cmd, stdout=PIPE, stderr=PIPE, shell=True)
        _, stderr = process.communicate()
        if stderr:
            print(stderr.decode("ascii"))
    

    def create_index(self) -> list[str]:
        self.target_dir = os.path.join(self.task_dir, "target/")
        targetpath = os.path.join(self.target_dir, "target.fasta" )
        os.makedirs(self.target_dir, exist_ok=True)
        if not os.path.exists(targetpath):
            if not os.path.exists("target.fasta"):
                raise FileNotFoundError(
                    "target.fasta is missing and was not found in target directory"
                )
            shutil.move("target.fasta", targetpath)

        cmd = f"cd {self.target_dir} && bowtie-build -f target.fasta {BASENAME} && cd ../../.."
        self.run_command(cmd=cmd)
        Index = []
        host_dir = os.path.join(self.task_dir,"host")
        for file in os.listdir(host_dir):
            print("file",file)
            host_fasta = os.path.join(host_dir, file)
            if os.path.isfile(host_fasta):
                index = os.path.splitext(file)[0]
                host_index_file = os.path.join("bowtie_index",host_dir,f"{index}.1.ebwt")
                if os.path.exists(host_fasta) and not os.path.exists(host_index_file):
                    cmd = f"cd {host_dir} && bowtie-build -f {index}.fasta {index} && cd ../.."
                    self.run_command(cmd=cmd)
                Index.append(index)
        return Index 
        
    def run_bowtie(self, index:str, fasta_path:str, output_path:str) -> None:
        cmd = f"bowtie -x {index} -a -f {fasta_path} -v 3 > {output_path}"
        self.run_command(cmd=cmd)
