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
text = ironbark

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
use_variable_symbol_output = True
move_to_front = True


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
    save_bit = 0 if save else 1
    print(f"--> 0{save_bit}: Output literal, {save_text[save_bit]}")
    to_write_len = len(to_write)
    len_bits = fib_encode(to_write_len)
    print(f"  > {len_bits} ({len(len_bits)} bits): length of {to_write_len} characters to write")
    print(f"  > (Output of {repr(''.join(to_write))} in binary): {to_write_len*8} bits")
    output_lines.append(
        f"0{save_bit} {len_bits} {repr(''.join(to_write))}"
    )
    to_write = []
    return (2 + len(len_bits) + to_write_len * 8)


def output_symbols_equal_length(output=True):
    """
    Output each symbol, assuming all symbols are of length equal to the
    number of bits needed to represent the largest symbol index.
    """
    # The symbols in to_write here are the actual strings that matched.
    global to_write
    output_bit = 0 if output else 1
    print(f"--> 1{output_bit}: Output symbols, {symbol_text[output_bit]}")
    to_write_len = len(to_write)
    len_bits = fib_encode(to_write_len)
    print(f"  > {len_bits} ({len(len_bits)} bits): length of {to_write_len} characters to write")
    symbol_bits = int(log2(len(symbols))) + 1
    encoded_list = []
    for symbol in to_write:
        encoded_list.append(f"{repr(symbol)}({symbols[symbol]['num']})")
        if move_to_front:
            move_symbol_to_front(symbol)
    symb_encoding = ','.join(encoded_list)
    print(f"  > Symbols {symb_encoding} in {symbol_bits} bits each = {to_write_len * symbol_bits} bits")
    for symbol in to_write:
        symbols[symbol]['used in output'] += 1
        # Increment all 2+ length sub-symbols as prefixes used in this symbol
        for l in range(2, len(symbol)):
            symbols[symbol[:l]]['used as prefix'] += 1
    output_lines.append(
        f"1{output_bit} {len_bits} {symb_encoding} in {symbol_bits} bits each"
    )
    to_write = []
    return (2 + len(len_bits) + to_write_len * symbol_bits)


def output_symbols_variable_length(output=True):
    """
    Output each symbol with a variable number of bits.
    """
    # The symbols in to_write here are the actual strings that matched.
    global to_write
    output_bit = 0 if output else 1
    print(f"--> 1{output_bit}: Output symbols, {symbol_text[output_bit]}")
    bit_count = 2
    to_write_len = len(to_write)
    len_bits = fib_encode(to_write_len)
    bit_count += len(len_bits)
    print(f"  > {len_bits} ({len(len_bits)} bits): {to_write_len} symbols to write")
    encoded_symbols = list()
    for symbol in to_write:
        symbol_d = symbols[symbol]
        sym_fib = fib_encode(symbol_d['num']+1)  # Make sure symbol 0 encoded!
        print(f"  > symbol {repr(symbol)}({symbol_d['num']}) => {sym_fib} ({len(sym_fib)} bits)")
        encoded_symbols.append(sym_fib)
        bit_count += len(sym_fib)
        symbol_d['used in output'] += 1
        # Increment all 2+ length sub-symbols as prefixes used in this symbol
        for l in range(2, len(symbol)):
            symbols[symbol[:l]]['used as prefix'] += 1
        if move_to_front:
            move_symbol_to_front(symbol)
    output_lines.append(
        f"1{output_bit} {len_bits} {' '.join(encoded_symbols)}"
    )
    to_write = []
    return bit_count


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
            if len(symbols) == 2<<int(log2(len(symbols))):
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
