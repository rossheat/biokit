import '../constants/strings.dart';
import '../helpers/validateNucleotideSequence.dart';

Map<String, dynamic> matchMotif(
    String nucleotideSequence, String motifSequence, String sequenceType) {
  Map<String, String> sequenceValidationMap =
      validateNucleotideSequence(nucleotideSequence, sequenceType);
  String validNucleotideSequence = sequenceValidationMap[kSequence];

  Map<String, String> motifValidationMap =
      validateNucleotideSequence(motifSequence, sequenceType);
  String validMotifSequence = motifValidationMap[kSequence];

  int motifSequenceLength = validMotifSequence.length;
  int nucleotideSequenceLength = validNucleotideSequence.length;

  if (motifSequenceLength > nucleotideSequenceLength) {
    throw ("Invalid Motif Sequence Length Error. The motif sequence must not be longer than than the primary sequence. Motif sequence length: $motifSequenceLength nt. Primary sequence length: $nucleotideSequenceLength nt.");
  }

  List<Map<String, int>> matchIndices = [];
  Map<String, dynamic> matchMotifResultMap = {};

  int matchCount = 0;
  for (var i = 0; i <= nucleotideSequenceLength - motifSequenceLength; i++) {
    if (validNucleotideSequence.substring(i, i + motifSequenceLength) ==
        validMotifSequence) {
      matchCount++;
      matchIndices
          .add({kStartIndex: i, kEndIndex: i + motifSequenceLength - 1});
    }
  }

  matchMotifResultMap[kMatchCount] = matchCount;
  matchMotifResultMap[kMatchIndices] = matchIndices;
  return matchMotifResultMap;
}
