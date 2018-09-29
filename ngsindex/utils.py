"""Supporting classes.
"""
from pathlib import Path
from struct import Struct
from typing import Sequence, Tuple, Union, Optional, Any, cast

from xphyle import xopen


STRUCT_TYPES = dict(
    char=("c", 1),
    byte=("b", 1),
    ubyte=("B", 1),
    bool=("?", 1),
    short=("h", 2),
    ushort=("H", 2),
    int=("i", 4),
    uint=("I", 4),
    long=("l", 4),
    ulong=("l", 4),
    ll=("q", 8),
    ull=("Q", 8),
    float=("f", 4),
    double=("d", 8),
)
"""Data types that can be read from a binary file. A dict that maps the data type
name to a tuple of (`format_char`, `num_bits`), where `format_char` is a valid format
character in the python `struct` module.
"""


class DefaultDict(dict):
    """Similar to a `collections.defaultdict` except that, in addition to
    returning a default value for missing keys, the key=default pair is also
    stored to the dict.

    Args:
        default: A default value, a function that returns the default value,
            or None.
        pass_key_as_arg: Whether to pass the key to the `default` function.
    """

    def __init__(
        self,
        *args,
        default: Optional[Any] = None,
        pass_key_as_arg: bool = False,
        **kwargs,
    ) -> None:
        super().__init__(*args, **kwargs)
        if not callable(default):
            self.default = lambda: default
        else:
            self.default = default
        self.pass_key_as_arg = pass_key_as_arg

    def __getitem__(self, key):
        if key in self:
            return super().__getitem__(key)
        else:
            if self.pass_key_as_arg:
                value = self.default(key)
            else:
                value = self.default()
            super().__setitem__(key, value)
            return value


class BinReader:
    """Read binary information from a file and covert it into python data
    types. Provides two methods each for most data types described in
    https://docs.python.org/3/library/struct.html:
    * read_<dtype>() - read a single value of the specified type and returns
      the value.
    * read_<dtype>s(n) - reads n values of the specified type and returns a
      list.

    `read_string(strlen)` reads a string (type 's') of `strlen` characters and
    returns it as a string type rather than an array of characters.

    `read_vector(types)` and `read_vectors(types, n)` are lower-level methods
    that enable reading of arbitrary formats.

    Either `path` or `fileobj` must be specified.

    Args:
        path: The file to read from.
        byte_order: The endian-ness; '<' for little-endian, '>' for big-endian.
        fileobj: The file object to read from.
    """

    def __init__(
        self, path: Optional[Path] = None, byte_order: str = "@", fileobj=None
    ) -> None:
        if not (path or fileobj):
            raise ValueError("'path' or 'fileobj' must be provided")
        self.path = path
        self._reader = None
        self._close = True
        self.byte_order = byte_order
        self._struct_cache = DefaultDict(default=Struct, pass_key_as_arg=True)

        if fileobj:
            self._reader = fileobj
            self._close = False
            if not self.path and hasattr(fileobj, "name"):
                self.path = fileobj.name

    def __enter__(self):
        if not self.is_open:
            self.open()
        return self

    @property
    def is_open(self) -> bool:
        """Whether the reader is currently open."""
        return self._reader is not None

    def open(self) -> None:
        """Opens the reader.
        """
        self._reader = xopen(self.path, "rb", context_wrapper=True)

    def __exit__(self, exception_type, exception_value, traceback):
        self.close()

    def close(self) -> None:
        """Closes the reader.
        """
        if self._reader:
            if self._close:
                self._reader.close()
            self._reader = None

    def __getattr__(self, name: str):
        if not name.startswith("read_"):
            raise AttributeError(name)
        name = name[5:]
        if name.endswith("s"):
            return lambda n: self._read_array(name[:-1], n=n)
        else:
            return lambda: self._read_array(name)[0]

    def read_string(
        self, strlen: int, encoding: Optional[Union[bool, str]] = None
    ) -> Union[str, bytes]:
        """Reads `strlen` bytes and returns either a byte or character string.

        Args:
            strlen: Number of bytes to read.
            encoding: Byte encoding, or True to use default encoding.

        Returns:
            If `encoding` is None the raw bytes are returned, otherwise the
            bytes are decoded and a string is returned.
        """
        return self.decode(self._read_array("s", n=strlen)[0], encoding)

    def read_strings(
        self,
        size: Optional[Union[int, Sequence[int]]] = None,
        n: Optional[int] = None,
        encoding: Optional[str] = None,
    ) -> Sequence[Union[str, bytes]]:
        """Convenience method for reading multiple strings.

        Args:
            size: String size, or sequence of string sizes.
            n: Number of strings to read. If None, number of strings is
                inferred from length of `size`.
            encoding: Byte encoding, or True to use default encoding.

        Returns:
            A list of strings.
        """
        if n:
            size_list = [size] * n
        elif isinstance(size, int):
            size_list = [size]
        else:
            size_list = cast(Sequence[int], size)
        total = sum(size_list)
        concat = self.decode(self._read_array("s", n=total)[0], encoding)
        start = 0
        strings = []
        for size in size_list:
            end = start + size
            strings.append(concat[start:end])
            start = end
        return strings

    @staticmethod
    def decode(bytestr: bytes, encoding: Union[bool, str] = True) -> Union[str, bytes]:
        """Decode `bytes`.

        Args:
            bytestr: The bytes to decode.
            encoding: A valid byte endcoding, or True to use default encoding.

        Returns:
            The decoded bytes, or the bytes themselves if `encoding` is None/False.
        """
        if encoding is True:
            return bytestr.decode()
        elif encoding:
            return bytestr.decode(encoding)
        else:
            return bytestr

    def _read_array(self, type_: str, n: int = 1) -> Sequence:
        """Reads an array of values (i.e. all of the same data type) of length
        `n` with total size `size` (in bytes).

        Args:
            type_: The data type (see python `struct` library).
            n: The number of values in the array.

        Returns:
            A tuple of values.
        """
        fmt, size = STRUCT_TYPES.get(type_, (type_, 1))
        return self.read_struct("{}{}{}".format(self.byte_order, n, fmt), size * n)

    def read_vector(self, types: Sequence[Tuple[str, int]]) -> Sequence:
        """Read multiple consecutive values, which may be of different types,
        and return them in a tuple.

        Args:
            types: Sequence of tuples, where each tuple is (type, n).

        Returns:
            A tuple of values the same length as `types`.
        """
        return self.read_vectors(types, n=1)[0]

    def read_vectors(self, types: Sequence[Tuple[str, int]], n: int) -> Sequence:
        """Read `n` vectors.

        Args:
            types: The types of the values in each vector.
            n: The number of values in the array.

        Returns:
            A list of tuples.
        """
        formats = []
        nbytes = 0
        for type_name, type_n in types:
            fmt, size = STRUCT_TYPES.get(type_name, (type_name, 1))
            formats.append("{}{}{}".format(self.byte_order, type_n, fmt))
            nbytes += size * type_n
        return self.read_structs("".join(formats), nbytes, n)

    def read_struct(self, fmt: str, nbytes: int) -> Sequence:
        """Reads a struct, which is defined by a format string using the
        mini-language described in the python `struct` package.

        Args:
            fmt: The format string.
            nbytes: The total size in bytes of the struct.

        Returns:
            A tuple of values.
        """
        return self.read_structs(fmt, nbytes, 1)[0]

    def read_structs(self, fmt: str, nbytes: int, n: int) -> Sequence:
        """Reads `n` structs of the given format, each having size of `nbytes`.

        Args:
            fmt: The format string.
            nbytes: The size of each struct.
            n: The number of structs to read.

        Returns:
            A list of tuples.

        Notes:
            This can be faster than calling read_struct() `n` times, as it
            reads all of the required bytes (nbytes * n) in a single read
            operation.
        """
        if not self.is_open:
            raise RuntimeError(
                "Must call 'open' on BinReader before calling any read methods."
            )
        total = nbytes * n
        bytestr = self._reader.read(total)
        if len(bytestr) < total:
            raise EOFError()
        return [
            self._struct_cache[fmt].unpack(bytestr[i:(i + nbytes)])
            for i in range(0, total, nbytes)
        ]
