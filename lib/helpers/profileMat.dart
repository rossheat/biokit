import 'package:biokit/constants/lists.dart';
import 'package:biokit/constants/strings.dart';
import 'package:biokit/helpers/validateEqLens.dart';
import 'package:biokit/helpers/validateNucSeq.dart';
import 'package:meta/meta.dart';

Map<String, List<int>> profileMat(
    {@required List<String> nucSeqs, @required String seqType}) {
  for (var nucSeq in nucSeqs) {
    validateNucSeq(nucSeq: nucSeq, seqType: seqType);
  }

  int nucSeqsLen = validateEqLens(nucSeqs: nucSeqs);

  Map<String, List<int>> nucFreq = {};
  for (var nuc in seqType == kDNA ? dnaNucs : rnaNucs) {
    for (var i = 0; i < nucSeqsLen; i++) {
      int nucCount = 0;
      for (var nucSeq in nucSeqs) {
        if (nucSeq[i] == nuc) {
          nucCount++;
        }
      }
      nucFreq[nuc] == null
          ? nucFreq[nuc] = [nucCount]
          : nucFreq[nuc].add(nucCount);
      nucCount = 0;
    }
  }

  return nucFreq;
}
