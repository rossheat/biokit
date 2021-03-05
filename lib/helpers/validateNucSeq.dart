import 'package:meta/meta.dart';
import '../constants/lists.dart';
import '../constants/strings.dart';
import '../errors/invalidSeqErrMsg.dart';
import 'validateSeqType.dart';

Map<String, String> validateNucSeq(
    {@required String nucSeq, @required String seqType}) {
  if (nucSeq.length < 3) {
    throw ('Invalid Sequence Length Error. Sequence has less than three elements.');
  }

  String upperNucSeq = nucSeq.toUpperCase();
  String lowerSeqType = validateSeqType(seqType: seqType);

  upperNucSeq.split('').asMap().forEach((idx, nuc) {
    if (lowerSeqType == kDNA
        ? !dnaNucs.contains(nuc)
        : !rnaNucs.contains(nuc)) {
      throw (invalidSeqErrMsg(monomer: nuc, index: idx, seqType: seqType));
    }
  });

  return {kSeq: upperNucSeq, kSeqType: lowerSeqType};
}
