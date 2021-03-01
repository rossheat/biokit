import '../constants/maps.dart';
import '../constants/strings.dart';
import '../helpers/reverseNucleotideSequence.dart';
import '../helpers/validateNucleotideSequence.dart';

String complementaryStrand(
    String nucleotideSequence, String sequenceType, bool reverse) {
  Map<String, String> validationMap = validateNucleotideSequence(
      nucleotideSequence = nucleotideSequence, sequenceType = sequenceType);

  String validNucleotideSequence = validationMap[kSequence];
  String validSequenceType = validationMap[kSequenceType];

  String complementarySequence = validNucleotideSequence
      .split('')
      .map((nucleotide) => validSequenceType == kDNA
          ? DNAComplementaryNucleotides[nucleotide]
          : RNAComplementaryNucleotides[nucleotide])
      .join();

  return reverse
      ? reverseNucleotideSequence(complementarySequence)
      : complementarySequence;
}
