import '../constants/strings.dart';
import '../helpers/validateNucleotideSequence.dart';

String DNAToRNA(String nucleotideSequence) {
  // Fix sequenceType = 'DNA'
  Map<String, String> validationMap =
      validateNucleotideSequence(nucleotideSequence = nucleotideSequence, kDNA);

  String validDNASequence = validationMap[kSequence];

  return validDNASequence.replaceAll(RegExp(r'T'), 'U');
}
