import 'package:meta/meta.dart';

String invalidSeqTypeErrMsg({@required String seqType}) {
  return "Invalid Sequence Type Error. '$seqType' is not a valid sequence type. Please enter the argument 'dna' or 'rna' for parameter 'seqType'.";
}
