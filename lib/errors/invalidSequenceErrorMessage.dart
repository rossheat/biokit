String invalidSequenceErrorMessage(
    String nucleotide, int index, String sequenceType) {
  return "Invalid ${sequenceType.toUpperCase()} Sequence Error. Character '$nucleotide' found at index position $index (zero-based) is not a valid ${sequenceType.toUpperCase()} nucleotide.";
}
