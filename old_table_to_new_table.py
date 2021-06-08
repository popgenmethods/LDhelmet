from __future__ import print_function
import argparse
import numpy
import sys

#version number (uint64_t) + num_conf_list (uint64_t) + theta (double)
NUM_HEADER_BITS = 64 + 64 + 64 

OLD_VERSION_BITS = numpy.array([0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,1,1,1,1,1,1,0,1,1,1,1,1])
NEW_VERSION_BITS = numpy.array([0,0,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,1,1,1,1,1,1,0,1,1,1,1,1])


if __name__=="__main__":
	parser = argparse.ArgumentParser(description='Converts tables created using LDhelmet versions 1.8 or lower to be compatible with LDhelmet v1.9+')
	parser.add_argument('-i', '--input', type=str, help='Filename of old LDhelmet style table',  required=True)
	parser.add_argument('-o', '--output', type=str, help='Filename for output',  required=True)
	parser.add_argument('--interpolate', help='Interpolate between values -- do not use if using an LDhelmet generated table.  LDpop or LDhat tables that had been converted to LDhelmet format are okay.', action='store_true')
	args = parser.parse_args()

	#get the bytes for whether or not to interpolate
	if args.interpolate:
		interpolate = numpy.unpackbits(numpy.ones((1), dtype="uint8"))
	else:
		interpolate = numpy.unpackbits(numpy.zeros((1), dtype="uint8"))

	#Load the bytes, convert to bits
	input_bytes = numpy.fromfile(args.input, dtype = "uint8")
	input_bits = numpy.unpackbits(input_bytes)

	#check versioning
	version_bits = input_bits[0:64]
	if not numpy.all(version_bits == OLD_VERSION_BITS):
		print("This converter only works on old tables.\n")
		sys.exit(0)
	input_bits[0:64] = NEW_VERSION_BITS

	#Add in the interpolate byte, write to file
	output_bits = numpy.hstack([input_bits[:NUM_HEADER_BITS], interpolate, input_bits[NUM_HEADER_BITS:]])
	output_bytes = numpy.packbits(output_bits)
	output_bytes.tofile(args.output)
	
