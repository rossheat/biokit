import 'package:biokit/constants/maps.dart';
import 'package:biokit/constants/strings.dart';
import 'package:biokit/helpers/validateNucSeq.dart';
import 'package:meta/meta.dart';

Map<String, int> codonFreq(
    {@required String nucSeq,
    @required String seqType,
    @required String aminoAcid}) {
  Map<String, String> validationMap =
      validateNucSeq(nucSeq: nucSeq, seqType: seqType);

  String validSeq = validationMap[kSeq];

  if (validSeq.length % 3 != 0) {
    throw ("Invalid {seqType.toUpperCase()} Sequence Length Error. The {seqType} sequence contains ${validSeq.length} nt. The {seqType.toUpperCase()} Sequence must contain a number of nt. divisible by three in order to be translated into an amino acid sequence");
  }

  Map<String, int> codonFreqMap = {};

  for (var i = 0; i < validSeq.length - 2; i += 3) {
    String codon = validSeq.substring(i, i + 3);
    String fetchedAminoAcid =
        seqType == kDNA ? dnaCodonToAA[codon] : rnaCodonToAA[codon];
    if (fetchedAminoAcid == aminoAcid.toUpperCase()) {
      codonFreqMap[codon] == null
          ? codonFreqMap[codon] = 1
          : codonFreqMap[codon]++;
    }
  }
  return codonFreqMap;
}
