// TODO - Add start codon functionality
// TODO - Report nucleotides before start codon and after stop codon

import '../constants/lists.dart';
import '../constants/maps.dart';
import '../constants/strings.dart';
import '../helpers/validateNucleotideSequence.dart';

Map<String, dynamic> RNAtoProtein(String nucleotideSequence) {
  Map<String, String> validationMap =
      validateNucleotideSequence(nucleotideSequence, kRNA);
  String validRNASequence = validationMap[kSequence];

  if (validRNASequence.length % 3 != 0) {
    // The sequence is not divisble by 3
    throw ("Invalid RNA Sequence Length Error. The RNA sequence contains ${validRNASequence.length} nt. The RNA Sequence must contain a number of nt. divisible by three in order to be translated into an amino acid sequence");
  }

  String polypeptideSequence = '';
  for (var i = 0; i < validRNASequence.length - 2; i += 3) {
    String codon = validRNASequence.substring(i, i + 3);

    if (translationStopCodons.contains(codon)) {
      return {
        kPolypeptideSequence: polypeptideSequence,
        kTranslationStopCodon: codon,
        kRNASequenceNucleotideCount: validRNASequence.length,
        kPolypeptideSequenceAminoAcidCount: polypeptideSequence.length
      };
    }
    polypeptideSequence +=
        RNACodonsToAminoAcids[validRNASequence.substring(i, i + 3)];
  }
  return {
    kPolypeptideSequence: polypeptideSequence,
    kTranslationStopCodon: null,
    kRNASequenceNucleotideCount: validRNASequence.length,
    kPolypeptideSequenceAminoAcidCount: polypeptideSequence.length
  };
}
