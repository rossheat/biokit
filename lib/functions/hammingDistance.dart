import '../constants/strings.dart';
import '../helpers/validateNucleotideSequence.dart';

int hammingDistance(String nucleotideSequenceOne, String nucleotideSequenceTwo,
    String sequenceType) {
  Map<String, String> validationMapOne =
      validateNucleotideSequence(nucleotideSequenceOne, sequenceType);
  Map<String, String> validationMapTwo =
      validateNucleotideSequence(nucleotideSequenceTwo, sequenceType);

  String validNucleotideSequenceOne = validationMapOne[kSequence];
  String validNucleotideSequenceTwo = validationMapTwo[kSequence];

  int sequenceOneLength = validNucleotideSequenceOne.length;
  int sequenceTwoLength = validNucleotideSequenceTwo.length;

  if (sequenceOneLength != sequenceTwoLength) {
    throw ('Unequal Sequence Length Error. Sequences must be of equal length to calculate Hamming Distance. The first sequence contains $sequenceOneLength nt., while the second sequence contains $sequenceTwoLength nt.');
  }

  int pointMutationCount = 0;

  validNucleotideSequenceOne
      .split('')
      .asMap()
      .forEach((index, sequenceOneNucleotide) {
    if (sequenceOneNucleotide != validNucleotideSequenceTwo[index]) {
      pointMutationCount++;
    }
  });
  return pointMutationCount;
}
