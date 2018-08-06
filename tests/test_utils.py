import random
import pytest
from ngsindex.utils import DefaultDict, BinReader


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
    with open(datadir / 'small.bam.bai', 'rb') as inp:
        br = BinReader(fileobj=inp, byte_order='<')
        with pytest.raises(AttributeError):
            br.foo()
    # TODO: test strings and vectors
