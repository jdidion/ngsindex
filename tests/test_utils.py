import io
import random
import struct

import pytest

from ngsindex.utils import BinReader, DefaultDict


def _generate_strings(n):
    for i in range(n):
        yield "".join(chr(random.randrange(32, 127)) for _ in range(10))


def _generate_key_value(n):
    for key in enumerate(_generate_strings(n)):
        value = random.random()
        yield (key, value)


def test_dict_ordered():
    d = {}
    keys = []
    for key, value in _generate_key_value(100):
        keys.append(key)
        d[key] = value
    assert keys == list(d.keys())


def test_DefaultDict():
    dd = DefaultDict(default=int)
    dd["a"] += 1
    assert dd["a"] == 1

    dd = DefaultDict(default=0)
    dd["a"] += 1
    assert dd["a"] == 1

    dd = DefaultDict(default=lambda k: k + 1, pass_key_as_arg=True)
    assert dd[1] == 2


def test_OrderedDefaultDict():
    odd = DefaultDict(default=int)
    keys = []
    values = []
    for key, value in _generate_key_value(100):
        keys.append(key)
        values.append(value)
        odd[key] += value
    assert keys == list(odd.keys())
    assert values == list(odd.values())


def test_BinReader(datadir):
    with pytest.raises(ValueError):
        BinReader()
    with open(datadir / "small.bam.bai", "rb") as inp:
        br = BinReader(fileobj=inp, byte_order="<")
        with pytest.raises(AttributeError):
            br.foo()


def test_BinReader_requires_open():
    # Reading before open() raises a RuntimeError rather than crashing.
    br = BinReader(fileobj=io.BytesIO(b"\x01\x00\x00\x00"), byte_order="<")
    br.close()
    with pytest.raises(RuntimeError):
        br.read_int()


def test_BinReader_name_from_fileobj():
    buf = io.BytesIO(b"\x00")
    buf.name = "from-fileobj.bin"
    br = BinReader(fileobj=buf)
    assert br.path == "from-fileobj.bin"


def test_BinReader_scalar_and_array():
    data = struct.pack("<i", 7) + struct.pack("<3i", 1, 2, 3)
    with BinReader(fileobj=io.BytesIO(data), byte_order="<") as br:
        assert br.read_int() == 7
        assert list(br.read_ints(3)) == [1, 2, 3]


def test_BinReader_eof():
    with BinReader(fileobj=io.BytesIO(b"\x01\x00"), byte_order="<") as br:
        with pytest.raises(EOFError):
            br.read_int()  # needs 4 bytes, only 2 available


def test_BinReader_read_string():
    with BinReader(fileobj=io.BytesIO(b"hello"), byte_order="<") as br:
        assert br.read_string(5, encoding=True) == "hello"
    with BinReader(fileobj=io.BytesIO(b"hello"), byte_order="<") as br:
        assert br.read_string(5) == b"hello"
    with BinReader(fileobj=io.BytesIO("héllo".encode()), byte_order="<") as br:
        assert br.read_string(6, encoding="utf-8") == "héllo"


def test_BinReader_read_strings():
    payload = b"AABBBCCCC"
    # Explicit sizes inferred from the sequence length.
    with BinReader(fileobj=io.BytesIO(payload), byte_order="<") as br:
        assert br.read_strings(size=[2, 3, 4], encoding=True) == ["AA", "BBB", "CCCC"]
    # n copies of a single size.
    with BinReader(fileobj=io.BytesIO(b"ABCD"), byte_order="<") as br:
        assert br.read_strings(size=2, n=2, encoding=True) == ["AB", "CD"]
    # Single int size, no n.
    with BinReader(fileobj=io.BytesIO(b"XY"), byte_order="<") as br:
        assert br.read_strings(size=2, encoding=True) == ["XY"]


def test_BinReader_vector_single_type():
    # read_vector / read_vectors compose a struct format from one (type, count)
    # entry -- this is how Bin.read_chunks consumes pairs of ull values.
    data = struct.pack("<QQ", 10, 20) + struct.pack("<QQ", 30, 40)
    with BinReader(fileobj=io.BytesIO(data), byte_order="<") as br:
        assert br.read_vector([("ull", 2)]) == (10, 20)
        assert br.read_vector([("ull", 2)]) == (30, 40)


def test_BinReader_read_vectors_n():
    data = struct.pack("<QQ", 1, 2) + struct.pack("<QQ", 3, 4)
    with BinReader(fileobj=io.BytesIO(data), byte_order="<") as br:
        assert br.read_vectors([("ull", 2)], n=2) == [(1, 2), (3, 4)]


def test_BinReader_read_struct_raw_type():
    # A type name not in STRUCT_TYPES falls through to a raw struct format char.
    data = struct.pack("<H", 513)
    with BinReader(fileobj=io.BytesIO(data), byte_order="<") as br:
        assert br.read_struct("<H", 2) == (513,)


def test_BinReader_seek_and_skip():
    data = struct.pack("<4i", 100, 200, 300, 400)
    with BinReader(fileobj=io.BytesIO(data), byte_order="<") as br:
        br.skip(4)  # skip first int
        assert br.read_int() == 200
        br.seek(0)
        assert br.read_int() == 100
        br.seek(12)  # jump to the 4th int
        assert br.read_int() == 400


def test_BinReader_context_manager_reopen():
    buf = io.BytesIO(b"\x05\x00\x00\x00")
    br = BinReader(fileobj=buf, byte_order="<")
    assert br.is_open
    with br as same:
        assert same is br
        assert br.read_int() == 5
