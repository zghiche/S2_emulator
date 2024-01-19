import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Numeric Input Parser")
    parser.add_argument("number", type=int, help="Enter a numeric value")
    return parser.parse_args()

def compress_value(value, exponent_bits=4, mantissa_bits=3, truncation_bits=0, lsb=1):
    value = value/lsb
    saturation_code = (1 << (exponent_bits + mantissa_bits)) - 1
    saturation_value = ((1 << (mantissa_bits + truncation_bits + 1)) - 1) << ((1 << exponent_bits) - 2)

    if value > saturation_value:
        return saturation_value

    bitlen = 0
    shifted_value = value >> truncation_bits
    valcopy = shifted_value
    while valcopy != 0:
        valcopy >>= 1
        bitlen += 1

    if bitlen <= mantissa_bits:
        compressed_code = shifted_value
        compressed_value = shifted_value << truncation_bits
        return compressed_value

    # Build exponent and mantissa
    exponent = bitlen - mantissa_bits
    mantissa = (shifted_value >> (exponent - 1)) & ~(1 << mantissa_bits)

    print(bin(exponent), bin(mantissa))
    print(mantissa<<exponent)
    compressed_value = ((1 << mantissa_bits) | mantissa) << (exponent - 1); 
    return compressed_value

if __name__ == "__main__":
    args = parse_args()
    value = args.number
    result = compress_value(value)
    
    print("Original value", value, bin(value))
    print("Result", result, bin(result))
