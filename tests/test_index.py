"""Test cases ported from biogo:
https://github.com/biogo/hts/blob/master/bam/bam_test.go
"""
import pytest
import gzip
import io

from . import IndexTestData, ChunkTestData, ChunkTestDataSet
from ngsindex import (
    Offset, Bin, Chunk, CsiRefIndex, BaiRefIndex, DualRefIndex,
    CoordinateIndex, BaiIndex, CsiIndex, TbiIndex, MissingIndexError, IndexType,
    resolve_index_file, parse_index
)


def test_resolve(datadir):
    assert resolve_index_file(datadir / 'small.bam', IndexType.BAI) == \
        datadir / 'small.bam.bai'
    assert resolve_index_file(datadir / 'small.vcf.gz', IndexType.CSI) == \
        datadir / 'small.vcf.gz.csi'
    assert resolve_index_file(datadir / 'small.vcf.gz', IndexType.TBI) == \
        datadir / 'small.vcf.gz.tbi'
    assert resolve_index_file(datadir / 'small.bam', IndexType.TBI, error=False) is None
    with pytest.raises(MissingIndexError):
        resolve_index_file(datadir / 'small.bam', IndexType.TBI, error=True)


def test_parse(datadir):
    bai_file = resolve_index_file(datadir / 'small.bam', IndexType.BAI)
    assert isinstance(parse_index(bai_file), BaiIndex)
    csi_file = resolve_index_file(datadir / 'small.vcf.gz', IndexType.CSI)
    assert isinstance(parse_index(csi_file), CsiIndex)
    tbi_file = resolve_index_file(datadir / 'small.vcf.gz', IndexType.TBI)
    assert isinstance(parse_index(tbi_file), TbiIndex)
    with pytest.raises(ValueError):
        parse_index(bai_file, IndexType.CSI)


def test_offset():
    with pytest.raises(ValueError):
        Offset()

    offset = Offset(file_offset=98, block_offset=7)
    assert offset.offset == 6422535
    assert str(offset) == 'Offset(98, 7)'

    offset = Offset(6422535)
    assert offset.file_offset == 98
    assert offset.block_offset == 7
    assert str(offset) == 'Offset(98, 7)'


def test_bin():
    offset = Offset(6422535)
    chunk = Chunk(begin=offset, end=offset)
    _bin = Bin(bin_num=1, chunks=[chunk])
    assert _bin.num_chunks == 1
    assert _bin[0] == chunk
    assert str(_bin) == 'Bin(1): Chunks=[Chunk((98, 7), (98, 7))]'


def test_bam_ref_index():
    assert list(DualRefIndex().iter_chunks()) == []

    offset = Offset(6422535)
    chunk = Chunk(begin=offset, end=offset)
    _bin = Bin(bin_num=1, chunks=[chunk])
    ref_index = DualRefIndex.empty(reference_id=0)
    ref_index.bins = [_bin]
    ref_index.intervals = [offset]
    assert ref_index.num_bins == 1
    assert ref_index.num_intervals == 1
    assert list(ref_index.iter_chunks()) == [(_bin, chunk)]


def test_bam_index():
    index = BaiIndex()
    with pytest.raises(ValueError):
        _ = index.num_ref

    offset = Offset(6422535)
    chunk = Chunk(begin=offset, end=offset)
    _bin = Bin(bin_num=1, chunks=[chunk])
    ref_index = DualRefIndex(reference_id=0, bins=[_bin], offsets=[offset])
    index = BaiIndex(version=1, ref_indexes=[ref_index])
    assert len(index) == 1
    assert index[0] == ref_index


def test_parse_bai(test_bai_index_data: IndexTestData):
    _test_parse(test_bai_index_data)


def test_csi_ref_index():
    _ = CsiRefIndex.empty()


def test_parse_csi(test_csi_index_data: IndexTestData):
    _test_parse(test_csi_index_data)


def test_parse_tbi(test_tbi_index_data: IndexTestData):
    tbi = _test_parse(test_tbi_index_data)
    assert not tbi.is_sam
    assert tbi.is_vcf
    assert not tbi.is_zero_based
    assert tbi.meta_char == '#'


def _test_parse(index_data: IndexTestData):
    if index_data.err:
        with pytest.raises(index_data.err):
            _ = index_data.parse()
    else:
        index = index_data.parse()
        assert index == index_data.expect
        return index


def _test_chunk(index: CoordinateIndex, test_data: ChunkTestData):
    ref = index[test_data.ref_index]
    actual = list(ref.iter_chunks())[test_data.begin:test_data.end]
    if actual is None and test_data.expect is None:
        return
    assert actual == test_data.expect


def _test_chunk_set(test_data_set: ChunkTestDataSet):
    data_bytes = bytes(test_data_set.data)
    if test_data_set.compressed:
        data = io.BytesIO(gzip.decompress(data_bytes))
    else:
        data = io.BytesIO(data_bytes)
    index = parse_index(fileobj=data)
    for test_data in test_data_set.tests:
        if test_data.err:
            with pytest.raises(test_data.err):
                _test_chunk(index, test_data)
        else:
            _test_chunk(index, test_data)


# These test cases are designed to test the translation
# between genomic coordinates and index chunks. Since we
# don't implement that yet, skipping these tests.

# def test_chunks():
#    for test_data in TEST_CHUNK_DATA:
#        with pytest.subTest(test_data.name):
#            pytest._test_chunk_set(test_data)


KB = 2 ** 10


# Build a mock index
chrom1_offsets = [
    Offset(file_offset=1 * KB, block_offset=0),
    Offset(file_offset=1 * KB, block_offset=30 * KB),  # 0KB + 50KB; 30 KB
    Offset(file_offset=1 * KB, block_offset=50 * KB),  # 0KB + 20KB; 20 KB
    Offset(file_offset=9 * KB, block_offset=10 * KB),  # 12KB - 40KB; 64-40 KB
    Offset(file_offset=21 * KB, block_offset=10 * KB),  # 8KB +  0KB; 64+0 KB
]


chrom2_offsets = []


chrom3_offsets = [
    Offset(file_offset=0, block_offset=0),  # skipped
    Offset(file_offset=21 * KB, block_offset=60 * KB),  # 0KB + 50KB; 50 KB
    Offset(file_offset=41 * KB, block_offset=20 * KB),  # 20KB - 40KB; 80-40 KB
    Offset(file_offset=41 * KB, block_offset=20 * KB),  # 0KB +  0KB; 0 KB
    Offset(file_offset=65 * KB, block_offset=50 * KB),  # 24KB + 30KB; 96+30 KB
]


chrom1 = BaiRefIndex(offsets=chrom1_offsets)
chrom2 = BaiRefIndex(offsets=chrom2_offsets)
chrom3 = BaiRefIndex(offsets=chrom3_offsets)
good_index = BaiIndex(version=1, num_ref=2, ref_indexes=[chrom1, chrom2, chrom3])


bad_index = BaiIndex(
    version=1, num_ref=1, ref_indexes=[BaiRefIndex(offsets=[chrom1_offsets[0]])]
)


# big interval gets split into 4 pieces
piece_size = 16 * KB // 4
split_intervals = [
    ("chr3", (48 * KB) + (i * piece_size), (48 * KB) + ((i + 1) * piece_size))
    for i in range(4)
]


expected_merged = (
    [
        ("chr1", 0, 32 * KB),
        ("chr1", 32 * KB, 48 * KB),
        ("chr1", 48 * KB, 56 * KB),
        ("chr1", 56 * KB, 64 * KB),
        ("chr1", 64 * KB, 100000),
        ("chr3", 16 * KB, 48 * KB),
    ]
    + split_intervals
    + [("chr3", 64 * KB, 100000)]
)


expected_merged_16kb = [
    ("chr1", 0, 64 * KB),
    ("chr1", 64 * KB, 100000),
    ("chr3", 16 * KB, 48 * KB),
    ("chr3", 48 * KB, 100000),
]
