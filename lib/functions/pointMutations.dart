import 'package:meta/meta.dart';
import '../constants/strings.dart';
import '../helpers/validateNucSeq.dart';

int pointMutations(
    {@required String nucSeqOne,
    @required String nucSeqTwo,
    @required String seqType}) {
  Map<String, String> validationMapOne =
      validateNucSeq(nucSeq: nucSeqOne, seqType: seqType);
  Map<String, String> validationMapTwo =
      validateNucSeq(nucSeq: nucSeqTwo, seqType: seqType);

  String validNucSeqOne = validationMapOne[kSeq];
  String validNucSeqTwo = validationMapTwo[kSeq];

  int seqOneLen = validNucSeqOne.length;
  int seqTwoLen = validNucSeqTwo.length;

  if (seqOneLen != seqTwoLen) {
    throw ('Unequal Sequence Length Error. Sequences must be of equal length to calculate Hamming Distance. The first sequence contains $seqOneLen nt., while the second sequence contains $seqTwoLen nt.');
  }

  int pointMutationCount = 0;

  validNucSeqOne.split('').asMap().forEach((index, seqOneNuc) {
    if (seqOneNuc != validNucSeqTwo[index]) {
      pointMutationCount++;
    }
  });
  return pointMutationCount;
}
