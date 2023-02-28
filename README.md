# LZPW

LZPW is an improvement on the Lempel-Ziv-Welch (LZW) compression algorithm,
which itself is an improvement on the Lempel-Ziv 78 (LZ78) algorithm.

The `lzpw` library implements this compression algorithm.

## Algorithm

The main description of the LZPW algorithm is in the `Algorithm.md` file.

Lempel-Ziv-Welch starts with a symbol table of 256 entries, one for each byte
value.  It then adds symbols as they are seen and outputs only symbols of the
longest possible symbol that matches the sequence of input bytes.

The problems with LZW, as I see it, are three fold:

Firstly, it is inefficient in encoding data it hasn't seen before.  Each byte
becomes at least nine bits and, in the standard LZW algorithm, becomes twelve
bits, a size increase of 50%.

Secondly, it doesn't quickly adapt to patterns.  A repeat of a long string has
to start by outputting the first two bytes as a symbol, then the next two, and
so forth.  While these represent a 25% reduction (16 bits represented in 12),
it would be good to find a way to encode this more efficiently.

We do want to remain true to the 'symbol' philosophy of LZ78 and LZW.
Therefore, the obvious solution to the second problem - to have a 'symbol'
that encodes 'repeat this many bytes from this far away in the raw data' - is
out as that is basically the LZ77 algorithm.  And so the solution to that is
to try and encode more recently used symbols in fewer bits.

The third problem is that the symbol table is bounded, usually at around 4095
symbols (in the original LZW fixed 12-bit symbol implementation).  The answer
to the symbol table filling is simply to flush all >1 byte symbols and start
again.  While this does allow the algorithm to 'adapt' to changes in symbol
frequency, it is expensive in loss of existing information.  It would be good
to find a way to 'prune' the symbol table as it grows.

So our goals are:

* Encode literal bytes efficiently.
* Encode symbols more efficiently.
* Keep the symbol table from growing too large.

And our corollary goal is to not re-implement LZ77.

# Concepts

## Variable length whole number encoding

One basic concept used in this algorithm is that of the []Fibonacci, or
Zeckendorf, representation](https://en.wikipedia.org/wiki/Fibonacci_coding)
of positive integers.

Basically, base radices are the Fibonacci numbers (1, 2, 3, 5, 8, 13, ...)
and numbers are written in little-endian bit order, with a one bit appended.
Since (due to Zeckendorf's Theorem) such a number cannot have two consecutive
one bits, reading two one bits is taken to signal the end of the number.  The
first ten positive integers ('Natural numbers') with their Fibonacci codings
are:

- 1 = 11
- 2 = 011
- 3 = 0011
- 4 = 1011
- 5 = 00011
- 6 = 10011
- 7 = 01011
- 8 = 000011
- 9 = 100011
- 10 = 010011

Lengths and symbol indices are encoded using this method.

## Symbols

The LZPW algorithm, like LZW and LZ78 before it, thinks in terms of 'symbols',
which are taken to be groups of two or more bytes.  In LZW and LZ78, symbols
are encoded using the same 'unit' of output as 'literal' bytes.  In LZPW,
symbols are differentiated from literal bytes, and can therefore be given
shorter encodings.

Each symbol is given a number, starting with 1 (because this is the smallest
number than can be encoded using the simplest Fibonacci coding).

## Move to front

Each time a symbol is used, it is moved to first place (i.e. given the number
1), and each symbol numbered between 1 and the used symbol's number is moved
up one place.  For example, if symbol number 5 is used, it becomes symbol 1;
symbol 1 becomes symbol 2, symbol 2 becomes symbol 3, symbol 3 becomes symbol
4, and symbol 4 becomes symbol 5.

By doing this, symbols that are used more frequently have lower symbol numbers
and therefore shorter encodings due to Fibonacci coding.

Another idea here is to not only move the recently used symbol to the front,
but all of the symbols which use it as a prefix, as well.  A symbol that uses
another as a prefix is an integral part of LZ77 and LZW encoding, and it
recognises that in the future we hope to use this longer symbol rather than
its shorter prefix.  Giving the 'suffix' symbols of this symbol smaller
numbers tries to reward the reuse of both prefix and suffix.

For example, if the symbol table has symbols for 'AB' and 'ABC', then when
'AB' is moved to number 1, 'ABC' is moved to number 2.  If 'ABC' was
previously symbol number 7, all symbols from 2 to 6 are now renumbered 3 to 7.

## Requirements

The `lzpw` library uses `pipenv` to install its requirements.

Python development libraries (providing `Python.h`) need to be installed for
the `bitstream` library.
