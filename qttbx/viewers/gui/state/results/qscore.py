from dataclasses import dataclass
from typing import List


from .results import Result

@dataclass(frozen=True)
class QscoreResult(Result):
  qscore_per_atom: List[float]

