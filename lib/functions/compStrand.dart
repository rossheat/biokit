import 'package:meta/meta.dart';
import '../constants/maps.dart';
import '../constants/strings.dart';
import '../helpers/reverseNucSeq.dart';
import '../helpers/validateNucSeq.dart';

String compStrand(
    {@required String nucSeq,
    @required String seqType,
    @required bool reversed}) {
  Map<String, String> validationMap =
      validateNucSeq(nucSeq: nucSeq, seqType: seqType);

  String validNucSeq = validationMap[kSeq];
  String validSeqType = validationMap[kSeqType];

  String compSeq = validNucSeq
      .split('')
      .map((nuc) => validSeqType == kDNA ? dnaCompNucs[nuc] : rnaCompNucs[nuc])
      .join();

  return reversed ? reverseNucSeq(nucSeq: compSeq) : compSeq;
}
