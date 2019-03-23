"""Parsers for various NGS indexing formats (tabix, BAI, CSI).

Notes:
    This is a direct port of the parsers in biogo:
    https://github.com/biogo/hts/blob/master/tabix/tabix.go
"""
from abc import ABCMeta, abstractmethod
from enum import Enum
from pathlib import Path
from typing import Dict, Iterator, Optional, Sequence, Tuple, Type, Union, cast

from ngsindex.utils import BinReader

# pylint: disable=protected-access
# noinspection PyProtectedMember
from ngsindex._version import get_versions

__version__ = get_versions()["version"]
del get_versions


METADATA_BIN_NUM = 37450
"""Magic number of bins in BAM indexes that contain chromosome metadata."""
TERMINATOR = b"\x00"
"""Terminator character used in Index file."""


Coordinates = Tuple[int, int]


class MissingIndexError(Exception):
    """The index file could not be resolved.
    """

    pass


class Offset:
    """A virtual file offset. Must specify either `offset` or both `file_offset` and
    `block_offset`.

    Args:
        offset: A 64 bit int, with the high-order 48 bits being the file offset and the
            low-order 16 bits being the block offset.
        file_offset: The offset into the file.
        block_offset: The offset into the uncompressed block that starts at that file
            offset.
    """

    def __init__(
        self,
        offset: Optional[int] = None,
        file_offset: Optional[int] = None,
        block_offset: Optional[int] = None,
    ) -> None:
        if offset is not None:
            self.offset = offset
            self.file_offset = offset >> 16
            self.block_offset = offset & 0xffff
        elif file_offset is not None and block_offset is not None:
            self.file_offset = file_offset
            self.block_offset = block_offset
            self.offset = (file_offset << 16) + block_offset
        else:
            raise ValueError("Must specify either offset or file_offset+block_offset")

    def as_tuple(self) -> Coordinates:
        """Gets the offset as a tuple.

        Returns:
            A tuple (file_offset, block_offset)
        """
        return self.file_offset, self.block_offset

    def __eq__(self, other: "Offset") -> bool:
        return self.offset == other.offset

    def __repr__(self) -> str:
        return f"Offset{self.as_tuple()}"


class Chunk:
    """A section of a Bin.

    Args:
        begin: The bin start.
        end: The bin end.
    """

    def __init__(self, begin: Union[Offset, int], end: Union[Offset, int]) -> None:
        if isinstance(begin, int):
            begin = Offset(begin)
        if isinstance(end, int):
            end = Offset(end)
        self.begin = begin
        self.end = end

    def as_tuple(self) -> Tuple[Coordinates, Coordinates]:
        """Gets the chunk as a tuple.

        Returns:
            A tuple (begin, end)
        """
        return self.begin.as_tuple(), self.end.as_tuple()

    def __eq__(self, other: "Chunk") -> bool:
        return self.as_tuple() == other.as_tuple()

    def __repr__(self) -> str:
        return f"Chunk{self.as_tuple()}"


class Bin:
    """A bin in the hierarchical index.

    Args:
        reader: The BinReader. If not None, this Bin is populated by reading
            the bin_num and chunks from the reader.
        bin_num: The sequential index of this bin.
        chunks: A sequence of chunks.
    """

    def __init__(
        self,
        reader: Optional[BinReader] = None,
        bin_num: Optional[int] = None,
        chunks: Optional[Sequence[Chunk]] = None,
    ) -> None:
        self.bin_num = bin_num
        self.chunks = chunks
        if reader:
            self.read_from(reader)

    def read_from(self, reader: BinReader) -> None:
        """Initializes this Bin from `reader`.

        Args:
            reader: The BinReader from which to read.
        """
        self.bin_num = reader.read_uint()
        num_chunks = reader.read_int()
        self.chunks = Bin.read_chunks(reader, num_chunks)

    @property
    def num_chunks(self) -> int:
        """The number of Chunks in this Bin.
        """
        return len(self.chunks) if self.chunks else 0

    def __getitem__(self, index: int) -> Chunk:
        return self.chunks[index]

    @staticmethod
    def read_chunks(reader: BinReader, num_chunks: int) -> Sequence[Chunk]:
        """Reads chunks from a BinReader. Can be overridden.

        Args:
            reader: The BinReader from which to read.
            num_chunks: The number of chunks to read.

        Returns:
            A sequence of Chunk objects.
        """
        return [
            Chunk(*chunk) for chunk in reader.read_vectors([("ull", 2)], n=num_chunks)
        ]

    def __eq__(self, other: "Bin") -> bool:
        return self.bin_num == other.bin_num and self.chunks == other.chunks

    def __repr__(self) -> str:
        return f"Bin({self.bin_num}): Chunks={self.chunks}"


class CsiBin(Bin):
    """A Bin in a CsiRefIndex.

    Args:
        reader: The BinReader from which to read.
        bin_num: The index of this bin.
        first_record_offset: The offset of the first record that overlaps this bin.
        chunks: The Chunks within this bin.
    """

    def __init__(
        self,
        reader: Optional[BinReader] = None,
        bin_num: Optional[int] = None,
        first_record_offset: Optional[int] = None,
        chunks: Optional[Sequence[Chunk]] = None,
    ):
        self.first_record_offset = first_record_offset
        super().__init__(reader, bin_num, chunks)

    def read_from(self, reader: BinReader) -> None:
        """Initializes this CsiBin from `reader`.

        Args:
            reader: The BinReader from which to read.
        """

        self.bin_num = reader.read_uint()
        self.first_record_offset = reader.read_ull()
        num_chunks = reader.read_int()
        self.chunks = CsiBin.read_chunks(reader, num_chunks)

    def __eq__(self, other: "CsiBin") -> bool:
        return (
            super().__eq__(other)
            and self.first_record_offset == other.first_record_offset
        )

    def __repr__(self) -> str:
        return f"CsiBin({self.bin_num}): Offset: {self.first_record_offset}, " \
               f"Chunks={self.chunks}"


class RefIndex(metaclass=ABCMeta):
    """Base for classes that store indexes for a single reference sequence.

    Args:
        reader: The BinReader from which to read.
        reference_id: The ID of the reference sequence.
        bins: The bins in the index.
    """

    def __init__(
        self,
        reader: Optional[BinReader] = None,
        reference_id: Optional[int] = None,
        bins: Optional[Sequence[Bin]] = None,
    ) -> None:
        self.reference_id = reference_id
        self.bins = bins
        if reader:
            self.read_from(reader)

    @classmethod
    def empty(cls, reference_id: Optional[int] = None) -> "RefIndex":
        """Creates an empty RefIndex.

        Args:
            reference_id: The ID of the reference sequence.

        Returns:
            An instance of a RefIndex subclass.
        """
        return cls(reference_id=reference_id, bins=[])

    def read_from(self, reader: BinReader) -> None:
        """Initializes this RefIndex from `reader`.

        Args:
            reader: The BinReader from which to read.
        """
        num_bins = reader.read_int()
        self.bins = self.read_bins(reader, num_bins)

    @property
    def num_bins(self) -> int:
        """The number of Bins in this index.
        """
        return len(self.bins) if self.bins else 0

    @abstractmethod
    def read_bins(self, reader: BinReader, num_bins: int) -> Sequence[Bin]:
        """Reads Bins from `reader`.

        Args:
            reader: The BinReader from which to read.
            num_bins: The number of Bins to read.

        Returns:
            A sequence of Bins.
        """
        pass  # pragma: no-cover

    def iter_chunks(self) -> Iterator[Tuple[Bin, Chunk]]:
        """Iterates over all chunks of all bins.

        Yields:
            Tuples of (bin, chunk)
        """
        if not self.bins:
            return
        for _bin in self.bins:
            if _bin.num_chunks > 0:
                for chunk in _bin.chunks:
                    yield (_bin, chunk)

    def __eq__(self, other: "RefIndex") -> bool:
        return self.reference_id == other.reference_id and self.bins == other.bins

    def __repr__(self) -> str:
        type_name = type(self).__name__
        return f"{type_name}({self.reference_id}): Bins={self.bins}"


class CsiRefIndex(RefIndex):
    """RefIndex implementation for CSI index format.
    """

    def read_bins(self, reader: BinReader, num_bins: int) -> Sequence[Bin]:
        return [CsiBin(reader) for _ in range(num_bins)]


class DualRefIndex(RefIndex):
    """A RefIndex that has both a binning index and a linear index.

    Args:
        reader: The BinReader from which to read.
        reference_id: The ID of the reference sequence.
        bins: The bins in the binning index.
        offsets: The interval offsets in the linear index.
    """

    def __init__(
        self,
        reader: Optional[BinReader] = None,
        reference_id: Optional[int] = None,
        bins: Optional[Sequence[Bin]] = None,
        offsets: Optional[Sequence[Offset]] = None,
    ) -> None:
        self.intervals = offsets
        super().__init__(reader, reference_id, bins)

    @classmethod
    def empty(cls, reference_id: Optional[int] = None) -> "DualRefIndex":
        """Creates an empty DualRefIndex.

        Args:
            reference_id: The ID of the reference sequence.

        Returns:
            An DualRefIndex object.
        """
        return cls(reference_id=reference_id, bins=[], offsets=[])

    @staticmethod
    def read_intervals(reader: BinReader, num_intervals: int) -> Sequence[Offset]:
        """Reads intervals from `reader`.

        Args:
            reader: The BinReader from which to read.
            num_intervals: The number of Offsets to read.

        Returns:
            A sequence of Offsets.
        """
        if num_intervals == 0:
            return []
        return [Offset(offset) for offset in reader.read_ulls(num_intervals)]

    @property
    def num_intervals(self) -> int:
        """The number of intervals in the index.
        """
        return len(self.intervals) if self.intervals else 0

    def read_from(self, reader: BinReader) -> None:
        """Initializes this DualRefIndex from `reader`.

        Args:
            reader: The BinReader from which to read.
        """
        super().read_from(reader)
        num_intervals = reader.read_int()
        self.intervals = DualRefIndex.read_intervals(reader, num_intervals)

    def read_bins(self, reader: BinReader, num_bins: int) -> Sequence[Bin]:
        return [Bin(reader) for _ in range(num_bins)]

    def __eq__(self, other: "DualRefIndex") -> bool:
        return super().__eq__(other) and self.intervals == other.intervals

    def __repr__(self) -> str:
        type_name = type(self).__name__
        return f"{type_name}({self.reference_id}): Bins={self.bins}, " \
               f"Intervals={self.intervals}"


class BaiRefIndex(DualRefIndex):
    """A DualRefIndex subclass for the BAI index format.

    Args:
        reader: The BinReader from which to read.
        reference_id: The ID of the reference sequence.
        unmapped_begin: Beginning offset of unmapped reads.
        unmapped_end: Ending offset of unmapped reads.
        num_mapped: Number of mapped reads.
        num_unmapped: Number of unmapped reads.
        kwargs: Additional keyword arguments to pass to the RefIndex constructor.
    """

    def __init__(
        self,
        reader: Optional[BinReader] = None,
        reference_id: Optional[int] = None,
        unmapped_begin: Optional[Coordinates] = None,
        unmapped_end: Optional[Coordinates] = None,
        num_mapped: Optional[int] = None,
        num_unmapped: Optional[int] = None,
        **kwargs,
    ) -> None:
        self.unmapped_begin = unmapped_begin
        self.unmapped_end = unmapped_end
        self.num_mapped = num_mapped
        self.num_unmapped = num_unmapped
        super().__init__(reader, reference_id, **kwargs)

    def read_bins(self, reader: BinReader, num_bins: int) -> Sequence[Bin]:
        bins = []
        for i in range(num_bins):
            b = Bin(reader)
            if b.bin_num == METADATA_BIN_NUM:
                self.unmapped_begin, self.unmapped_end = b.chunks[0].as_tuple()
                self.num_mapped, self.num_unmapped = b.chunks[1].as_tuple()
            else:
                bins.append(b)
        return bins


class Index(metaclass=ABCMeta):
    """Base class for index formats.

    Args:
        reader: The BinReader from which to read.
        version: Version of the index specification.
    """

    def __init__(
        self,
        reader: Optional[BinReader] = None,
        version: Optional[int] = None
    ) -> None:
        self.version = version
        if reader:
            self.read_from(reader)

    @classmethod
    def index_type(cls) -> "IndexType":
        return IndexType.from_class(cls)

    def read_from(self, reader: BinReader) -> None:
        """Initializes this Index from `reader`.

        Args:
            reader: The BinReader from which to read.
        """
        self.version = reader.read_byte()

    @abstractmethod
    def summarize(self) -> dict:
        """Summarize this index and return a dict that can be serialized to JSON.
        """
        pass

    def __eq__(self, other: "Index") -> bool:
        return (
            self.index_type() == other.index_type()
            and self.version == other.version
        )

    def __repr__(self) -> str:
        type_name = self.index_type().name
        return f"{type_name}(v{self.version})"


class CoordinateIndex(Index):
    """Base class for formats that index by coordinate.

    Args:
        reader: The BinReader from which to read.
        version: Version of the index specification.
        num_ref: Number of reference sequences.
        ref_indexes: Sequence of RefIndex objects.
        num_unmapped: Number of unmapped reads.
    """

    def __init__(
        self,
        reader: Optional[BinReader] = None,
        version: Optional[int] = None,
        num_ref: Optional[int] = None,
        ref_indexes: Optional[Sequence[RefIndex]] = None,
        num_unmapped: Optional[int] = None,
    ) -> None:
        self._num_ref = num_ref
        self.ref_indexes = ref_indexes
        self.num_unmapped = num_unmapped
        super().__init__(reader, version)

    @property
    def num_ref(self) -> int:
        """The number of reference sequences.
        """
        if self._num_ref is None:
            if self.ref_indexes is not None:
                self._num_ref = len(self.ref_indexes)
            else:
                raise ValueError("num_ref is unset")
        return self._num_ref

    def __len__(self) -> int:
        return self.num_ref

    def read_from(self, reader: BinReader) -> None:
        """Initializes this Index from `reader`.

        Args:
            reader: The BinReader from which to read.
        """

        super().read_from(reader)

        # Each type of index has its own header format.
        self._num_ref = self.parse_header(reader)

        self.ref_indexes = [
            self.create_ref_index(reader, i) for i in range(self.num_ref)
        ]

        # Unmapped reads is optional
        self.num_unmapped = None
        try:
            self.num_unmapped = reader.read_ull()
        except EOFError:
            pass

    def __getitem__(self, idx: int) -> RefIndex:
        return self.ref_indexes[idx]

    @abstractmethod
    def parse_header(self, reader: BinReader) -> int:
        """Parses the header portion of the index.

        Args:
            reader: The BinReader from which to read.

        Returns:
            Integer, number of reference sequences.
        """
        pass  # pragma: no-cover

    @abstractmethod
    def create_ref_index(self, reader: BinReader, reference_id: int) -> RefIndex:
        """Creates a RefIndex object and initalizes it from `reader`.

        Args:
            reader: The BinReader from which to read.
            reference_id: The ID of the reference sequence.

        Returns:
            A RefIndex object.
        """
        pass  # pragma: no-cover

    def summarize(self) -> dict:
        """Summarize this index and return a dict that can be serialized to JSON.
        """
        return {
            "num_references": self.num_ref,
            "num_unmapped": self.num_unmapped
        }

    def __eq__(self, other: "CoordinateIndex") -> bool:
        return (
            super().__eq__(other)
            and self.num_ref == other.num_ref
            and self.ref_indexes == other.ref_indexes
            and self.num_unmapped == other.num_unmapped
        )

    def __repr__(self) -> str:
        type_name = type(self).__name__
        return f"{type_name}(v{self.version}): {self.ref_indexes}"


class BaiIndex(CoordinateIndex):
    """Index subclass for BAI format.

    See:
        * https://samtools.github.io/hts-specs/SAMv1.pdf
    """

    def parse_header(self, reader: BinReader) -> int:
        # The only thing in the BAI header is the number of reference seqs
        return reader.read_int()

    def create_ref_index(self, reader: BinReader, reference_id: int) -> RefIndex:
        return BaiRefIndex(reader, reference_id)

    def summarize(self) -> dict:
        """Create a summary equivalent to samtools idxstats.

        Returns:
            Summary dict.
        """
        summary = super().summarize()
        summary["read_counts"] = []
        for ref_idx in self.ref_indexes:
            bai_ref_idx = cast(BaiRefIndex, ref_idx)
            summary["read_counts"].append(
                (bai_ref_idx.num_mapped, bai_ref_idx.num_unmapped)
            )
        return summary


class TbiIndex(CoordinateIndex):
    """Index subclass for TBI format.

    Args:
        file_format:
        col_seq:
        col_begin:
        col_end:
        meta:
        skip:
        ref_names:
        name_map:
        kwargs:

    See:
        * https://samtools.github.io/hts-specs/tabix.pdf
        * https://academic.oup.com/bioinformatics/article/27/5/718/262743

    TODO:
        Add reg2bin and reg2bins functions from the tabix spec.
    """

    def __init__(
        self,
        file_format: Optional[int] = None,
        col_seq: Optional[int] = None,
        col_begin: Optional[int] = None,
        col_end: Optional[int] = None,
        meta: Optional[int] = None,
        skip: Optional[int] = None,
        ref_names: Optional[Sequence[str]] = None,
        name_map: Optional[Dict[str, int]] = None,
        **kwargs
    ) -> None:
        self.file_format = file_format
        self.col_seq = col_seq
        self.col_begin = col_begin
        self.col_end = col_end
        self.meta = meta
        self.skip = skip
        self.ref_names = ref_names
        self.name_map = name_map
        super().__init__(**kwargs)

    def parse_header(self, reader: BinReader) -> int:
        # First 8 fields are int32s (see spec for details)
        (
            num_ref,
            self.file_format,
            self.col_seq,
            self.col_begin,
            self.col_end,
            self.meta,
            self.skip,
            len_name_concat,
        ) = reader.read_ints(8)

        if num_ref > 0:
            # Last header field is a string consisting of concatenated,
            # zero-terminated reference sequence names.
            ref_names_concat: bytes = reader.read_string(len_name_concat)

            if ref_names_concat[-1:] != TERMINATOR:
                raise ValueError(f"Invalid ref names: {ref_names_concat}")

            self.ref_names = tuple(
                ref.decode() for ref in ref_names_concat[:-1].split(TERMINATOR)
            )

            if num_ref != len(self.ref_names):
                raise ValueError(
                    f"num_seq field {num_ref} does not equal number of "
                    f"reference names {len(self.ref_names)}"
                )

            self.name_map = dict((n, i) for i, n in enumerate(self.ref_names))

        return num_ref

    def create_ref_index(self, reader: BinReader, reference_id: int) -> RefIndex:
        return DualRefIndex(reader, reference_id)

    @property
    def is_sam(self) -> bool:
        """Whether this is an index for a SAM/BAM file.
        """
        return self.file_format == 1

    @property
    def is_vcf(self) -> bool:
        """Whether this is an index for a VCF/BCF file.
        """
        return self.file_format == 2

    @property
    def is_zero_based(self) -> bool:
        """Whether the coordinates in this index are zero-based.
        """
        return (self.file_format & 0x10000) != 0

    @property
    def meta_char(self) -> str:
        """The meta character as a string.
        """
        return chr(self.meta)

    def __eq__(self, other: "TbiIndex") -> bool:
        return (
            super().__eq__(other)
            and self.file_format == other.file_format
            and self.col_seq == other.col_seq
            and self.col_begin == other.col_begin
            and self.col_end == other.col_end
            and self.meta == other.meta
            and self.skip == other.skip
            and self.ref_names == other.ref_names
            and self.name_map == other.name_map
        )


class CsiIndex(CoordinateIndex):
    """Index subclass for CSI format.

    Args:
        min_shift:
        depth:
        aux_len:
        aux:
        kwargs: Additional keyword arguments to pass to the Index constructor.

    See:
        * https://samtools.github.io/hts-specs/CSIv2.pdf

    TODO:
        Add reg2bin and reg2bins functions from the CSI spec.
    """

    def __init__(
        self,
        min_shift: Optional[int] = None,
        depth: Optional[int] = None,
        aux: Optional[Sequence[int]] = None,
        **kwargs
    ) -> None:
        self.min_shift = min_shift
        self.depth = depth
        self.aux = aux
        super().__init__(**kwargs)

    def parse_header(self, reader: BinReader) -> int:
        (self.min_shift, self.depth, aux_len) = reader.read_ints(3)

        self.aux = reader.read_ubytes(aux_len)

        return reader.read_int()

    def create_ref_index(self, reader: BinReader, reference_id: int) -> RefIndex:
        return CsiRefIndex(reader, reference_id)

    def __eq__(self, other: "CsiIndex") -> bool:
        return (
            super().__eq__(other)
            and self.min_shift == other.min_shift
            and self.depth == other.depth
            and self.aux == other.aux
        )


class SbiIndex(Index):
    """Index subclass for SBI format. Note that this is a proposed format based
    on Hadoop SBI - it is not yet accepted as part of the SAM specification.

    Args:
        file_length:
        md5:
        uuid:
        num_records:
        granularity:
        num_offsets:
        kwargs: Additional keyword arguments to pass to the Index constructor.

    See:
        * https://github.com/tomwhite/hts-specs/blob/sbi/SAMv1.tex
    """

    def __init__(
        self,
        file_length: Optional[int] = None,
        md5: Optional[bytes] = None,
        uuid: Optional[bytes] = None,
        num_records: Optional[int] = None,
        granularity: Optional[int] = None,
        num_offsets: Optional[int] = None,
        offsets: Optional[Sequence[Offset]] = None,
        **kwargs
    ):
        self.file_length = file_length
        self.md5 = md5
        self.uuid = uuid
        self.num_records = num_records
        self.granularity = granularity
        self.num_offsets = num_offsets
        self.offsets = offsets
        super().__init__(**kwargs)

    def read_from(self, reader: BinReader) -> None:
        super().read_from(reader)

        self.file_length = reader.read_ull()
        self.md5 = reader.read_string(16)
        self.uuid = reader.read_string(16)
        self.num_records = reader.read_ull()
        self.granularity = reader.read_ull()
        self.num_offsets = reader.read_ull()

        self.offsets = [
            Offset(offset)
            for offset in reader.read_ulls(self.num_offsets)
        ]

    def summarize(self) -> dict:
        return {
            "file_length": self.file_length,
            "num_records": self.num_records,
            "granularity": self.granularity,
            "num_offsets": self.num_offsets
        }


class IndexType(Enum):
    """Enumeration of supported index types.
    """

    BAI = "BAI", BaiIndex
    CSI = "CSI", CsiIndex
    TBI = "TBI", TbiIndex
    SBI = "SBI", SbiIndex

    @property
    def name(self) -> str:
        return self.value[0]

    @property
    def index_class(self) -> Type[CoordinateIndex]:
        return self.value[1]

    @staticmethod
    def from_name(name: str) -> "IndexType":
        for index_type in IndexType:
            if index_type.name == name:
                return index_type
        else:
            raise ValueError(f"No index type '{name}'")

    @staticmethod
    def from_class(cls) -> "IndexType":
        for index_type in IndexType:
            if index_type.index_class == cls:
                return index_type
        else:
            raise ValueError(f"No index type for class '{cls}'")


def resolve_and_parse_index(primary_file: Path, index_type: IndexType) -> Index:
    """Shortcut for
    parse_index(resolve_index_file(parimary_file, index_type), index_type)

    Args:
        primary_file: The primary file with the index to resolve.
        index_type: The type of index.

    Returns:
        An Index object.
    """
    index_file = resolve_index_file(primary_file, index_type, error=True)
    return parse_index(index_file, index_type)


def resolve_index_file(
    primary_file: Path, index_type: Optional[IndexType] = None,
    index_file: Optional[Path] = None, error: bool = False, create: bool = False
) -> Union[Path, None]:
    """Resolve the path of an index file that accompanies `primary_file`.

    Args:
        primary_file: The primary file to resolve.
        index_type: The type of index to resolve. If None, `index_file` must
            not be None.
        index_file: The index file. If None, the index file is assumed to be
            <primary_file>.<ext>.
        create: Whether to create the index if it does not exist.
        error: Whether to raise an error if the index file does not exist.

    Raises:
        MissingIndexError if the index file does not exist and `error` is
        True.
    """
    if not index_file:
        if index_type:
            index_file = Path(f"{str(primary_file)}.{index_type.name.lower()}")
        else:
            raise ValueError(
                "'index_type' must be specified if 'index_file' is None"
            )
    if not index_file.exists() and create:
        try:
            import pysam
            pysam.index(str(index_file))
        except:
            pass
    if not index_file.exists():
        if error:
            raise MissingIndexError(f"Index file {index_file} does not exist")
        else:
            return None
    return index_file


def parse_index(
    path: Optional[Path] = None, index_type: Optional[IndexType] = None, **kwargs
) -> Index:
    """Parse an index file and return an `Index` subclass.

    Args:
        path: Path to the index.
        index_type: The index type. TBI == generic tabix index,
            BAI == BAM index, CSI == coordinate-sorted index.
            This is only used to validate that the index file is
            of the expected type, since the true type is always
            detected from the file's magic string.
        kwargs: Keyword arguments to pass through to BinReader constructor.
    """
    with BinReader(path=path, byte_order="<", **kwargs) as reader:
        # First three bytes are a 'magic' value (the index type).
        magic = reader.read_string(3, True)
        try:
            detected_type = IndexType.from_name(magic)
        except ValueError:
            raise ValueError(f"Unsupported index type: {magic}")

        if index_type not in (None, detected_type):
            raise ValueError(
                f"Expected index type {index_type} does not match actual type {magic}"
            )

        return detected_type.index_class(reader=reader)
