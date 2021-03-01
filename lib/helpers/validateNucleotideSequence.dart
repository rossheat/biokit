import '../constants/lists.dart';
import '../constants/strings.dart';
import '../errors/invalidSequenceErrorMessage.dart';
import 'validateSequenceType.dart';

Map<String, String> validateNucleotideSequence(
    String nucleotideSequence, String sequenceType) {
  if (nucleotideSequence.length < 3) {
    throw ('Invalid Sequence Length Error. Sequence has less than three elements.');
  }

  String upperNucleotideSequence = nucleotideSequence.toUpperCase();

  String lowerSequenceType = validateSequenceType(sequenceType);

  upperNucleotideSequence.split('').asMap().forEach((index, nucleotide) {
    if (lowerSequenceType == kDNA
        ? !DNANucleotides.contains(nucleotide)
        : !RNANucleotides.contains(nucleotide)) {
      throw (invalidSequenceErrorMessage(
          nucleotide = nucleotide, index = index, sequenceType = sequenceType));
    }
  });

  return {kSequence: upperNucleotideSequence, kSequenceType: lowerSequenceType};
}
