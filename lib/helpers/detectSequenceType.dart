import '../constants/strings.dart';
import 'validateNucleotideSequence.dart';

Map<String, String> detectSequenceType(String nucleotideSequence) {
  bool isValidDNASequence;
  bool isValidRNASequence;
  String invalidDNAErrorMessage;
  String invalidRNAErrorMessage;

  try {
    validateNucleotideSequence(nucleotideSequence, kDNA);
    isValidDNASequence = true;
  } catch (errorMessage) {
    isValidDNASequence = false;
    invalidDNAErrorMessage = errorMessage;
  }

  try {
    validateNucleotideSequence(nucleotideSequence, kRNA);
    isValidRNASequence = true;
  } catch (errorMessage) {
    isValidRNASequence = false;
    invalidRNAErrorMessage = errorMessage;
  }

  if (isValidDNASequence == false && isValidRNASequence == false) {
    // Sequence is not DNA or RNA
    return {
      kResult: kInvalidSequence,
      kInvalidDNAErrorMessage: invalidDNAErrorMessage,
      kInvalidRNAErrorMessage: invalidRNAErrorMessage
    };
  } else if (isValidDNASequence == true && isValidRNASequence == false) {
    // Sequence is DNA, not RNA
    return {kResult: kDNA, kInvalidRNAErrorMessage: invalidRNAErrorMessage};
  } else if (isValidDNASequence == false && isValidRNASequence == true) {
    // Sequence is RNA, not DNA
    return {kResult: kRNA, kInvalidDNAErrorMessage: invalidDNAErrorMessage};
  } else {
    // Cannot tell if sequence is DNA or RNA (No T or U present in sequence)
    return {kResult: kAmbiguousSequence};
  }
}
