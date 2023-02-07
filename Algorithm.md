= The LZPW Algorithm =

The fundamental algorithm is in one of two states:

* Literal mode, where consecutive bytes are little-endian encoded.
* Symbol mode, where symbols from the symbol table are encoded.

Each symbol in the symbol table is two or more bytes that have been found
consecutively in the raw input.  Symbols can be prefixes for other symbols,
but the smallest-length symbol is two bytes and its prefix is a single byte.
The algorithm also remembers 'first bytes' (prefixes?).

The algorithm starts in literal mode.

The first byte is read, stored as the 'last sequence read', and added to the
list of literal bytes to be output.

Each subsequent byte read can be seen as a symbol when appended to the 'last
sequence read'.

If this symbol is not in the symbol table, the subsequent byte is added to
the list of literal bytes to be output.  The new symbol is stored and given
an index number.

== Adaptations ==

* Instead of bytes, assume a fixed bit width - e.g. colour palettes, DNA and
  protein encoding.
* Instead of bytes, assume variable bit width tokens (e.g. Unicode characters)
  read from an encoding, e.g. UTF-8, and treated as units that we call 'bytes'.
  This may already be done by the Python implementation.
