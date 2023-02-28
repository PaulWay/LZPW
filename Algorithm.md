# The LZPW Algorithm

The fundamental algorithm is in one of two states:

* Literal mode, where consecutive bytes are little-endian encoded.
* Symbol mode, where symbols from the symbol table are encoded.

The simplest way of describing the algorithm is: we write the longest literal
string that contains no groups that are in the symbol table.  Then we write
symbols for as long as we can group the input characters into symbols.  The
algorithm alternates between these two states.

Each symbol in the symbol table is two or more bytes that have been found
consecutively in the raw input.  Symbols can be prefixes for other symbols,
but the smallest-length symbol is two bytes and its prefix is a single byte.
The algorithm also remembers 'first bytes' (prefixes?).

The algorithm starts in literal mode.

The output stream consists of a sequence of blocks.  Each block starts with a
Fibonacci coded positive integer giving the length of this block.  Literal
blocks then list the bytes in the literal block, encoded as is.  Symbol blocks
then list each symbol's (current) index number, encoded in Fibonacci coding,
with the first symbol being given the index `1`.  Blocks alternate between
literal and then symbol blocks - blocks of the same type are always
concatenated together.

The first byte is read, stored as the 'last sequence read', and added to the
list of literal bytes to be output.

Each subsequent byte read can be seen as a potential symbol when appended to
the 'last sequence read'.

If this potential symbol is not in the symbol table, the new symbol is stored
and given the next highest index number.  In literal mode the subsequent byte
is added to the list of literal bytes to be output.  In symbol mode, we
remember the


## Possible variations

* Instead of bytes, assume a fixed range of symbols - e.g. colour palettes,
  DNA and protein encoding.  Doesn't have to be a power of two, although
  literal symbols would still have a fixed base-two encoding.
* Instead of bytes, assume variable bit width tokens (e.g. Unicode characters)
  read from an encoding, e.g. UTF-8, and treated as units that we call 'bytes'.
  This may already be done by the Python implementation.
