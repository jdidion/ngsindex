"""Additional tests targeting previously-uncovered branches in the public API.

These complement test_index.py (which is ported from biogo) by exercising the
SBI parser, the summarize() methods, the IndexType lookups, error paths in
resolve_index_file / parse_index, and the convenience helpers.
"""

import io
import json
import struct
from pathlib import Path

import pytest

from ngsindex import (
    BaiIndex,
    BaiRefIndex,
    Bin,
    Chunk,
    CsiIndex,
    DualRefIndex,
    IndexType,
    MissingIndexError,
    Offset,
    SbiIndex,
    TbiIndex,
    parse_index,
    resolve_and_parse_index,
    resolve_index_file,
)

SBI_OFFSETS = [6422535, 0, 100]


def _build_sbi(
    file_length=123456,
    md5=b"\x11" * 16,
    uuid=b"\x22" * 16,
    num_records=7,
    granularity=2,
    offsets=tuple(SBI_OFFSETS),
):
    """Build a spec-conformant SBI index blob.

    See https://github.com/tomwhite/hts-specs/blob/sbi/SAMv1.tex
    """
    blob = b"SBI"
    blob += struct.pack("<B", 1)  # version
    blob += struct.pack("<Q", file_length)
    blob += md5
    blob += uuid
    blob += struct.pack("<Q", num_records)
    blob += struct.pack("<Q", granularity)
    blob += struct.pack("<Q", len(offsets))
    for off in offsets:
        blob += struct.pack("<Q", off)
    return blob


SBI_DATA = _build_sbi()

# Real BAI/CSI/TBI fixtures live alongside test_index.py.
INDEX_DATA_DIR = Path(__file__).parent / "test_index"


def test_offset_empty():
    assert Offset(0).empty is True
    assert Offset(6422535).empty is False


def test_parse_sbi():
    index = parse_index(fileobj=io.BytesIO(SBI_DATA))
    assert isinstance(index, SbiIndex)
    assert index.version == 1
    assert index.file_length == 123456
    assert index.md5 == b"\x11" * 16
    assert index.uuid == b"\x22" * 16
    assert index.num_records == 7
    assert index.granularity == 2
    assert index.num_offsets == 3
    assert [o.offset for o in index.offsets] == SBI_OFFSETS


def test_sbi_summarize():
    index = parse_index(fileobj=io.BytesIO(SBI_DATA))
    summary = index.summarize()
    assert summary == {
        "file_length": 123456,
        "num_records": 7,
        "granularity": 2,
        "num_offsets": 3,
    }
    # Confirm it is JSON-serializable (this is what the CLI relies on).
    json.loads(json.dumps(summary))


def test_index_type_round_trip():
    for index_type in IndexType:
        assert IndexType.from_name(index_type.name) is index_type
        assert IndexType.from_class(index_type.index_class) is index_type


def test_index_type_from_name_unknown():
    with pytest.raises(ValueError):
        IndexType.from_name("ZZZ")


def test_index_type_from_class_unknown():
    with pytest.raises(ValueError):
        IndexType.from_class(int)


def test_index_index_type_helper():
    # Index.index_type() routes through IndexType.from_class.
    assert BaiIndex.index_type() is IndexType.BAI
    assert CsiIndex.index_type() is IndexType.CSI
    assert TbiIndex.index_type() is IndexType.TBI
    assert SbiIndex.index_type() is IndexType.SBI


def test_index_repr():
    assert repr(BaiIndex(version=3, num_ref=0, ref_indexes=[])) == "BaiIndex(v3): []"


def test_bai_iter_intervals():
    offset_a = Offset(file_offset=1, block_offset=0)
    offset_b = Offset(file_offset=2, block_offset=0)
    ref0 = DualRefIndex(reference_id=0, bins=[], offsets=[offset_a])
    ref1 = DualRefIndex(reference_id=1, bins=[], offsets=[offset_b])
    index = BaiIndex(version=1, num_ref=2, ref_indexes=[ref0, ref1])
    assert list(index.iter_intervals()) == [offset_a, offset_b]


def test_bai_summarize_read_counts():
    chunk = Chunk(begin=Offset(0), end=Offset(0))
    ref = BaiRefIndex(
        reference_id=0,
        num_mapped=8,
        num_unmapped=1,
        bins=[Bin(bin_num=1, chunks=[chunk])],
        offsets=[],
    )
    index = BaiIndex(version=1, num_ref=1, ref_indexes=[ref], num_unmapped=1)
    summary = index.summarize()
    assert summary["num_references"] == 1
    assert summary["num_unmapped"] == 1
    assert summary["read_counts"] == {0: [8, 1]}


def test_coordinate_index_getitem_and_len():
    ref = DualRefIndex(reference_id=0, bins=[], offsets=[])
    index = BaiIndex(version=1, num_ref=1, ref_indexes=[ref])
    assert len(index) == 1
    assert index[0] is ref


def test_resolve_index_file_requires_type_or_file(tmp_path):
    with pytest.raises(ValueError):
        resolve_index_file(tmp_path / "reads.bam")


def test_resolve_index_file_explicit_index_file(tmp_path):
    bam = tmp_path / "reads.bam"
    idx = tmp_path / "custom.bai"
    idx.write_bytes(b"")
    assert resolve_index_file(bam, index_file=idx) == idx


def test_resolve_index_file_missing_returns_none(tmp_path):
    assert resolve_index_file(tmp_path / "reads.bam", IndexType.BAI) is None


def test_resolve_index_file_missing_raises(tmp_path):
    with pytest.raises(MissingIndexError):
        resolve_index_file(tmp_path / "reads.bam", IndexType.BAI, error=True)


def test_resolve_index_file_create_best_effort(tmp_path):
    # `create=True` attempts pysam.index but swallows any failure; with no
    # pysam / no real BAM the file still does not exist, so None is returned.
    bam = tmp_path / "reads.bam"
    assert resolve_index_file(bam, IndexType.BAI, create=True) is None


def test_resolve_and_parse_index():
    index = resolve_and_parse_index(INDEX_DATA_DIR / "small.bam", IndexType.BAI)
    assert isinstance(index, BaiIndex)


def test_resolve_and_parse_index_missing(tmp_path):
    with pytest.raises(MissingIndexError):
        resolve_and_parse_index(tmp_path / "reads.bam", IndexType.TBI)


def test_parse_index_unsupported_magic():
    with pytest.raises(ValueError, match="Unsupported index type"):
        parse_index(fileobj=io.BytesIO(b"ZZZ\x00\x00\x00"))


def test_parse_index_type_mismatch():
    with pytest.raises(ValueError, match="does not match"):
        parse_index(fileobj=io.BytesIO(SBI_DATA), index_type=IndexType.BAI)


def test_parse_index_explicit_matching_type():
    index = parse_index(fileobj=io.BytesIO(SBI_DATA), index_type=IndexType.SBI)
    assert isinstance(index, SbiIndex)


def _build_tbi(num_ref, names_concat):
    """Build a minimal TBI blob with `num_ref` references and a header-only body.

    The trailing reference-index records are omitted; parse_header is reached
    and validates the names block before any ref index is read, which is all
    these error-path tests need.
    """
    blob = b"TBI"
    blob += struct.pack("<B", 1)  # version
    # num_ref, file_format, col_seq, col_begin, col_end, meta, skip, len_names
    blob += struct.pack("<8i", num_ref, 2, 1, 2, 0, 35, 0, len(names_concat))
    blob += names_concat
    return blob


def test_tbi_ref_names_missing_terminator():
    # Names block that does not end in a NUL terminator is rejected.
    data = _build_tbi(1, b"19")
    with pytest.raises(ValueError, match="Invalid ref names"):
        parse_index(fileobj=io.BytesIO(data))


def test_tbi_ref_names_count_mismatch():
    # num_ref claims 3 but only two names are present.
    data = _build_tbi(3, b"19\x0020\x00")
    with pytest.raises(ValueError, match="does not equal number"):
        parse_index(fileobj=io.BytesIO(data))


def test_repr_smoke():
    # The __repr__ methods build descriptive strings; exercise each so a
    # formatting typo in a debug repr is caught.
    offset = Offset(file_offset=98, block_offset=7)
    chunk = Chunk(begin=offset, end=offset)
    _bin = Bin(bin_num=1, chunks=[chunk])
    assert "Offset" in repr(offset)
    assert "Chunk" in repr(chunk)
    assert "Bin" in repr(_bin)
    dual = DualRefIndex(reference_id=0, bins=[_bin], offsets=[offset])
    assert "DualRefIndex(0)" in repr(dual)
    assert "Intervals" in repr(dual)
    bai = BaiIndex(version=1, num_ref=1, ref_indexes=[dual])
    assert "BaiIndex(v1)" in repr(bai)
    sbi = parse_index(fileobj=io.BytesIO(SBI_DATA))
    assert repr(sbi) == "SBI(v1)"


def test_index_equality_across_types():
    # Index.__eq__ compares index_type() and version.
    bai = BaiIndex(version=1, num_ref=0, ref_indexes=[])
    sbi = SbiIndex(version=1)
    assert bai != sbi
    other_bai = BaiIndex(version=2, num_ref=0, ref_indexes=[])
    assert bai != other_bai


def test_iter_chunks_skips_empty_bins():
    # A bin with no chunks must be skipped without yielding (branch 296->295).
    offset = Offset(6422535)
    chunk = Chunk(begin=offset, end=offset)
    full = Bin(bin_num=1, chunks=[chunk])
    empty = Bin(bin_num=2, chunks=[])
    ref = DualRefIndex(reference_id=0, bins=[empty, full, empty], offsets=[])
    assert list(ref.iter_chunks()) == [(full, chunk)]


def test_tbi_zero_references():
    # num_ref == 0 skips the names block entirely (branch 680->700).
    data = _build_tbi(0, b"")
    index = parse_index(fileobj=io.BytesIO(data))
    assert isinstance(index, TbiIndex)
    assert index.num_ref == 0
    assert index.ref_names is None
    assert index.name_map is None
