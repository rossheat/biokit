import 'package:biokit/helpers/validateNucSeq.dart';
import 'package:meta/meta.dart';
import '../constants/maps.dart';
import '../constants/strings.dart';

Map<String, dynamic> nucSeqToAA(
    {@required String nucSeq, @required String seqType, startIdx = 0}) {
  Map<String, String> validationMap =
      validateNucSeq(nucSeq: nucSeq, seqType: seqType);
  String validSeq = validationMap[kSeq];

  if (startIdx - 3 > validSeq.length) {
    throw ("Invalid Start Index Error. [startIdx] must be at least three less than the length of the sequence passed to [nucSeq]");
  }

  String aaSeq = '';
  for (var i = startIdx; i < validSeq.length - 2; i += 3) {
    String codon = validSeq.substring(i, i + 3);
    aaSeq += seqType == kDNA ? dnaCodonToAA[codon] : rnaCodonToAA[codon];
  }

  return {
    kAASeq: aaSeq,
    kSeqNucCount: validSeq.length - startIdx - 1,
    kAASeqCount: aaSeq.length
  };
}
