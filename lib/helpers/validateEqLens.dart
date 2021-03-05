import 'package:meta/meta.dart';

int validateEqLens({@required List<String> nucSeqs}) {
  if (nucSeqs.length < 2) {
    throw ('Invalid Number of Sequences Error. Must pass at least two sequences in order to compare lengths.');
  }
  int nucSeqLen = nucSeqs[0].length;
  for (var nucSeq in nucSeqs) {
    if (nucSeq.length != nucSeqLen) {
      throw ("Unequal Sequence Lengths Error. Sequences must be of equal length. Sequence '{nucSeqs[1]} is of length {nucSeqLen}, while is sequence {nucSeq} is of length {nucSeq.length}.");
    }
  }
  return nucSeqLen;
}
