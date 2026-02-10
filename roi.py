from dataclasses import dataclass, field



ROI_OFFSET = 400

@dataclass 
class ROI:
    name: str
    target_start: int
    target_end: int
    expected_structure_sequence: str

    roi_start: int = field(init=False)
    roi_end: int = field(init=False)
    roi_sequence: str = field(init=False) 

    def __post_init__(self) -> None:
        self.roi_start = max(self.target_start-ROI_OFFSET, 0)
        self.roi_end = min(self.target_end+ ROI_OFFSET, len(self.expected_structure_sequence))
        self.roi_sequence = self.expected_structure_sequence[self.roi_start: self.roi_end]

    def __repr__(self) -> str:
        buffer = min(ROI_OFFSET, self.target_start)
        seq_id = (f"SEQUENCE_ID={self.name}\n")
        seq_templ = f"SEQUENCE_TEMPLATE={self.roi_sequence}\n"
        seq_target = f"SEQUENCE_TARGET={buffer},{self.target_end-self.target_start}\n"
        return f"{seq_id}{seq_templ}{seq_target}=\n"