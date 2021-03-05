import 'package:meta/meta.dart';

String invalidSeqErrMsg(
    {@required String monomer, @required int index, @required String seqType}) {
  return "Invalid ${seqType.toUpperCase()} Sequence Error. Character '$monomer' found at index position $index (zero-based) is not a valid ${seqType.toUpperCase()} monomer.";
}
