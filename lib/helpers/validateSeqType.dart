import 'package:meta/meta.dart';
import '../constants/strings.dart';
import '../errors/invalidSeqTypeErrMsg.dart';

String validateSeqType({@required String seqType}) {
  String lowerSeqType = seqType.toLowerCase();

  if (lowerSeqType != kDNA && lowerSeqType != kRNA) {
    throw (invalidSeqTypeErrMsg(seqType: seqType));
  }

  return lowerSeqType;
}
