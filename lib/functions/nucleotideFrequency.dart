import '../constants/strings.dart';
import '../helpers/validateNucleotideSequence.dart';

Map<String, int> nucleotideFrequency(String nucleotideSequence, sequenceType) {
  Map<String, String> validationMap = validateNucleotideSequence(
      nucleotideSequence = nucleotideSequence, sequenceType = sequenceType);

  String validNucleotideSequence = validationMap[kSequence];

  Map<String, int> nucleotideFrequencyMap = {};
  validNucleotideSequence.split('').forEach((nucleotide) {
    if (nucleotideFrequencyMap.containsKey(nucleotide)) {
      nucleotideFrequencyMap[nucleotide]++;
    } else {
      nucleotideFrequencyMap[nucleotide] = 1;
    }
  });

  return nucleotideFrequencyMap;
}
