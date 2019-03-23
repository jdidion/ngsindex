from dataclasses import dataclass
import io
from typing import Callable, List, Any

from ngsindex import *


@dataclass
class IndexTestData:
    name: str
    data: Sequence[int]
    expect: CoordinateIndex
    err: Exception = None

    def parse(self):
        data = io.BytesIO(bytes(self.data))
        return parse_index(fileobj=data)


@dataclass
class ChunkTestData:
    begin: int
    end: int
    expect: Sequence[Chunk] = None
    ref_index: int = 0
    err: Exception = None


@dataclass
class ChunkTestDataSet:
    name: str
    data: Sequence[int]
    tests: Sequence[ChunkTestData]
    compressed: bool = False


def sparse_list(size: int, factory: Callable, *args: Tuple[int, Any]) -> List:
    sl = [factory(i) for i in range(size)]
    for idx, value in args:
        sl[idx] = value
    return sl
