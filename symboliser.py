#!/usr/bin/env python

from math import log, log2, sqrt

# Constants
ironbark = """It was the man from Ironbark that struck the Sydney town,
He wandered over street and park, he wandered up and down.
He loitered here, he loitered there, 'til he was like to drop,
Until at last in sheer despair he sought a barber shop.
"'Ere, shave me beard and whiskers off, I'll be a man of mark,
I'll go and do the Sydney toff up home in Ironbark"."""
simple = 'abababababababababababaa'
four = 'abcdabcdabcdabcdabcdabcdabcdabcd'
wander = 'wander, wandering wanderer - wondering wonder wand, wander.'

# with open('symboliser.py') as fh:
#     text = ''.join(fh.read())
text = wander

phi = (1+sqrt(5))/2
save_text = ['save output to symbol table', 'discard output']
symbol_text = ['output symbols', 'delete symbols']

fibs = [1, 2, 3]
fib_encodings = {1: '11', 2: '011', 3: '0011'}

def fib_encode(i):
    assert i > 0, "Cannot encode zero in Fibonacci encoding"
    if i in fib_encodings:
        return fib_encodings[i]
    # Find first fib not less than i; calculate and store if necessary
    fib_gt_idx = 2  # Don't start from end, i might be within array
    while fibs[fib_gt_idx] <= i:
        fib_gt_idx +=1
        if fib_gt_idx == len(fibs):
            fibs.append(fibs[-1] + fibs[-2])
    # Now calculate back from most significant bit
    bits = '1'  # to terminate
    for fib_idx in reversed(range(fib_gt_idx)):  # count from largest index down...
        if i >= fibs[fib_idx]:
            bits = '1' + bits
            i -= fibs[fib_idx]
        else:
            bits = '0' + bits
    fib_encodings[i] = bits
    return bits


# State information
symbols = {}
first_chars = set()
last_symbol = ''
to_write = []
state = 'literal'
bits_output = 0
mabulate_flag = False
output_lines = []

# Different processing ideas that we're trying out
use_variable_symbol_output = True
move_to_front = True
extra_output_bit = False
block_type_bit = False


def new_symbol(symbol):
    symbols[symbol] = {
        'num': len(symbols),
        'used as prefix': 0,
        'used in output': 0,
    }
    first_chars.add(symbol[0])


def move_symbol_to_front(symbol):
    """
    Move the symbol given to the front, incrementing the symbol number of
    all symbols less than this symbol's number.  We have to do this each time
    a symbol is written, because it then changes the symbol numbering and
    a subsequent symbol in the list of output symbols may have its number
    changed.

    Unfortunately at the moment, with only a dictionary holding the symbols
    and no other way to iterate through them, we have to do a full dictionary
    scan for this.
    """
    # Does it make a difference that we're choosing to output this list
    # after all the symbols have been emitted?  Could we output 'ab', 'cd'
    # and then 'cd', at which time 'cd' would be symbol 0 and have encoding
    # 1?
    symbol_d = symbols[symbol]
    source_num = symbol_d['num']
    print(f"    Moving symbol {repr(symbol)} from {source_num} to 0")
    if symbol_d['num'] == 0:
        return
    # moved_symbols = []
    for s_key, s_val in symbols.items():
        if s_val['num'] < source_num:
            # print(f"   -> mv {repr(s_key)} from {s_val['num']} to {s_val['num']+1}")
            # moved_symbols.append(repr(s_key))
            s_val['num'] += 1
    # print(f"    moved symbols {','.join(moved_symbols)}")
    symbol_d['num'] = 0


def output_literals(save=True):
    global to_write
    block_start = ''
    block_start_desc = ''
    if block_type_bit:
        block_start = '0'
        if extra_output_bit:
            save_bit = 0 if save else 1
            block_start += f"{save_bit}"
            block_start_desc = ', ' + save_text[save_bit]
    if block_start:
        block_start += ' '
    print(f"--> Output literals: {block_start}{block_start_desc}")
    bits_written = len(block_start)
    to_write_len = len(to_write)
    len_bits = fib_encode(to_write_len)
    bits_written += len(len_bits)
    print(f"  > {len_bits} ({len(len_bits)} bits): length of {to_write_len} characters to write")
    output_bits = " ".join(f"{ord(b):08b}" for b in to_write)
    output_bits_len = len(output_bits.replace(' ', ''))
    bits_written += output_bits_len
    print(f"  > {output_bits} {repr(''.join(to_write))}: {output_bits_len} bits")
    output_lines.append(
        f"l {block_start}{len_bits} {output_bits} {repr(''.join(to_write))}"
    )
    to_write = []
    print(f"  = {bits_written} bits")
    return bits_written


def output_symbols_equal_length(output=True):
    """
    Output each symbol, assuming all symbols are of length equal to the
    number of bits needed to represent the largest symbol index.
    """
    # The symbols in to_write here are the actual strings that matched.
    global to_write

    block_start = ''
    block_start_desc = ''
    if block_type_bit:
        block_start = '0'
        if extra_output_bit:
            output_bit = 0 if output else 1
            block_start += f"{output_bit}"
            block_start_desc = ', ' + symbol_text[output_bit]
    if block_start:
        block_start += ' '
    print(f"--> Output symbols: {block_start}{block_start_desc}")
    bits_written = len(block_start)

    to_write_len = len(to_write)
    len_bits = fib_encode(to_write_len)
    bits_written += len_bits
    print(f"  > {len_bits} ({len(len_bits)} bits): length of {to_write_len} characters to write")
    symbol_bits = int(log2(len(symbols))) + 1
    encoded_symbols = list()
    encoded_reprs = list()
    for symbol in to_write:
        encoded_symbols.append(f"{symbols[symbol]['num']:0{symbol_bits}b}")
        encoded_reprs.append(symbol)
        if move_to_front:
            move_symbol_to_front(symbol)
    bits_written += to_write_len * symbol_bits
    symb_encoding = ' '.join(encoded_symbols)
    print(f"  > Symbols {symb_encoding} in {symbol_bits} bits each = {to_write_len * symbol_bits} bits")
    for symbol in to_write:
        symbols[symbol]['used in output'] += 1
        # Increment all 2+ length sub-symbols as prefixes used in this symbol
        for l in range(2, len(symbol)):
            symbols[symbol[:l]]['used as prefix'] += 1
    output_lines.append(
        f"s {block_start}{len_bits} {symb_encoding} {encoded_reprs}"
    )
    to_write = []
    print(f"  = {bits_written} bits")
    return bits_written


def output_symbols_variable_length(output=True):
    """
    Output each symbol with a variable number of bits.
    """
    # The symbols in to_write here are the actual strings that matched.
    global to_write
    block_start = ''
    block_start_desc = ''
    if block_type_bit:
        block_start = '0'
        if extra_output_bit:
            output_bit = 0 if output else 1
            block_start += f"{output_bit}"
            block_start_desc = ', ' + symbol_text[output_bit]
    if block_start:
        block_start += ' '
    print(f"--> Output symbols: {block_start}{block_start_desc}")
    bits_written = len(block_start)

    to_write_len = len(to_write)
    len_bits = fib_encode(to_write_len)
    bits_written += len(len_bits)
    print(f"  > {len_bits} ({len(len_bits)} bits): {to_write_len} symbols to write")
    encoded_symbols = list()
    encoded_reprs = list()
    for symbol in to_write:
        symbol_d = symbols[symbol]
        sym_fib = fib_encode(symbol_d['num']+1)  # Make sure symbol 0 encoded!
        print(f"  > symbol {repr(symbol)}({symbol_d['num']}) => {sym_fib} ({len(sym_fib)} bits)")
        encoded_symbols.append(sym_fib)
        encoded_reprs.append(symbol)
        bits_written += len(sym_fib)
        symbol_d['used in output'] += 1
        # Increment all 2+ length sub-symbols as prefixes used in this symbol
        for l in range(2, len(symbol)):
            symbols[symbol[:l]]['used as prefix'] += 1
        if move_to_front:
            move_symbol_to_front(symbol)
    output_lines.append(
        f"s {block_start}{len_bits} {' '.join(encoded_symbols)} {encoded_reprs}"
    )
    to_write = []
    print(f"  = {bits_written} bits")
    return bits_written


def find_all_sym_matches(current_symbols: list):
    """
    Yield lists of symbol sequences that match as much of the input list as
    we can.  This starts with the minimum length (2) and looks for increasing
    length symbols until it can find no more here.  We use this technique to
    recursively search for further matches.

    The complication here is that we get the list of currently matched
    symbols as a list of strings, rather than one entire string.  Rather than
    try and mangle symbols and step through them character by character, we
    combine it into a big string and then work through that.
    """
    symbol_string = ''.join(current_symbols)
    symbol_string_len = len(symbol_string)
    # Find all the (remaining) suffix lists that have symbols in our symbol
    # table from this start index
    def find_suffixes_from(start_idx: int):
        end_idx = start_idx + 2
        while end_idx <= symbol_string_len:
            this_sym = symbol_string[start_idx:end_idx]
            if this_sym not in symbols:
                break
            for suffixes in find_suffixes_from(end_idx):
                yield [this_sym] + suffixes
            yield [this_sym]
            end_idx += 1

    # Return the list of symbol lists, where the last symbol is right at the
    # end of the symbol string.  This is roughly in order of longest symbols
    # last, but pays no attention to the encoding length of each symbol.
    return [
        symbol_list
        for symbol_list in find_suffixes_from(0)
        if symbol_list[-1] == symbol_string[-len(symbol_list[-1]):]
    ]


def output_symbols(output=True):
    if use_variable_symbol_output:
        return output_symbols_variable_length(output)
    else:
        return output_symbols_equal_length(output)


for char in text:
    print(f"<-- {state}: last {repr(last_symbol)}+{repr(char)} ", end='')
    if not last_symbol:
        last_symbol = char
        to_write.append(char)
        print("first character")
        continue
    sought_symbol = last_symbol + char
    if sought_symbol not in symbols:
        if state == 'symbols' and to_write:
            print(f"??? symbol list: {sorted(symbols.keys())}")
            print(f"??? symbols matched: {to_write}, {repr(last_symbol)}, next char {repr(char)}")
            print(f"??? possible other symbol matches to get us to this point: {find_all_sym_matches(to_write + [last_symbol] + [char])}")
        new_symbol(sought_symbol)
        print(f"not a symbol, add symbol {symbols[sought_symbol]['num']:3d} (mabulate {mabulate_flag})")
        if state == 'symbols':
            # If we've read 'abcdab' and we've just read 'c', we don't have an
            # 'abc' symbol yet, but we might be able to write symbols 'ab' and
            # 'cd' if 'c' is the prefix to a symbol.
            if not mabulate_flag:
                print(f"... mabulate flag false, saving last matched symbol {repr(last_symbol)}")
                to_write.append(last_symbol)
            # If the symbol table is at a power of two size, then it may be
            # more efficient to output the current symbols at the current bit
            # length, and then the subsequent symbols can be at the new bit
            # length.  Need to work out what the efficient changeover length
            # is.
            if len(symbols) == 2<<int(log2(len(symbols))) and not use_variable_symbol_output:
                print(f"... Symbol table is an even power of two, maybe we should write out the symbols?")
            # We might have got 'ca', where 'c' is a prefix (and matched
            # last time we got here) and 'a' is a prefix but 'ca' is not yet
            # in the symbol table.  We'll know because the 'prefix' should be
            # the last character of the last symbol. I haven't found a good
            # name for this yet, but we need to hold the 'c' and look for the
            # next character as a prefix. I call it the mabulate flag.  If
            # that's false, we've got a prefix ('c') and we're looking for a
            # symbol.  If it's true, got the 'c', the 'a' is a prefix, but
            # this is part of the 'we didn't find a symbol' code so we can
            # say that this is not a symbol and we need to go back to
            # literals.
            if char in first_chars and not mabulate_flag:
                print(f"... {repr(char)} is a prefix, looking for next symbol, mabulate flag -> True")
                mabulate_flag = True
                last_symbol = char
            else:
                # Have to switch to a literal block
                bits_output += output_symbols()
                state = 'literal'
                if mabulate_flag:
                    print(f"... mabulate flag True, also add {repr(last_symbol)}")
                    to_write.append(last_symbol)
                mabulate_flag = False
                print(f"... mabulate flag -> False, go to literal mode and add {repr(char)}")
                to_write.append(char)
        else:
            if mabulate_flag:
                print("!!! mabulate flag True, in literal mode?  mabulate flag -> False")
                mabulate_flag = False
            print(f"... add character {repr(char)}")
            to_write.append(char)
        # Finally, this is the start of the next prefix we try to match.
        last_symbol = char
    else:
        mabulate_flag = False
        print(f"matches symbol {repr(sought_symbol)} ({symbols[sought_symbol]['num']}) (mabulate -> False)")
        # always more efficient to encode to symbols at this point
        if state == 'literal':
            # The last byte in the literals is the prefix of this symbol, so
            # we have to remove it from the output.
            print(f"... removing {repr(to_write[-1])} from literals: prefix of new symbol")
            to_write.pop()
            bits_output += output_literals()
            state = 'symbols'
        # We don't append the symbol here - because this might be the prefix
        # of a longer symbol.  We append the symbol when we no longer match,
        # and we're in symbol mode, and the mabulate flag isn't set (I think?)

        last_symbol = sought_symbol

print(f"End of text.  State={state}, mabulate_flag={mabulate_flag}, sought_symbol={repr(sought_symbol)}, last_symbol={repr(last_symbol)}, char={repr(char)}")
if state == 'symbols':
    if not mabulate_flag:
        to_write.append(sought_symbol)
    bits_output += output_symbols()

if mabulate_flag:
    print(f"Un-mabulate {repr(last_symbol)} for output.")
    to_write.append(last_symbol)
if to_write:
    bits_output += output_literals()

print(f"Total input bits: {len(text)*8}")
print(f"Total output bits: {bits_output}")
print(f"Compression ratio: {len(text)*8/bits_output:0.4f}:1")

print("Symbols:")
for symbol in sorted(symbols, key=lambda s: symbols[s]['num']):
    sym_data = symbols[symbol]
    print(
        f"{sym_data['num']:3d}: {repr(symbol)} "
        f"output {sym_data['used in output']} times, "
        f"prefixed {sym_data['used as prefix']} times."
    )

print("Output:")
for line in output_lines:
    print(line)
