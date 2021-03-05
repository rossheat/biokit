import 'package:biokit/constants/lists.dart';
import 'package:biokit/constants/strings.dart';
import 'package:biokit/errors/invalidSeqErrMsg.dart';
import 'package:meta/meta.dart';

String validateAASeq({@required String aaSeq}) {
  String upperAASeq = aaSeq.toUpperCase();

  upperAASeq.split('').asMap().forEach((idx, aa) {
    if (!aminoAcids.contains(aa)) {
      throw (invalidSeqErrMsg(monomer: aa, index: idx, seqType: kAA));
    }
  });

  return upperAASeq;
}
