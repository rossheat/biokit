import 'package:meta/meta.dart';
import '../constants/strings.dart';
import '../helpers/validateNucSeq.dart';

Map<String, int> nucFreq({@required String nucSeq, @required String seqType}) {
  Map<String, String> validationMap =
      validateNucSeq(nucSeq: nucSeq, seqType: seqType);

  String validNucSeq = validationMap[kSeq];

  Map<String, int> nucFreqMap = {};
  validNucSeq.split('').forEach((nuc) {
    if (nucFreqMap.containsKey(nuc)) {
      nucFreqMap[nuc]++;
    } else {
      nucFreqMap[nuc] = 1;
    }
  });

  return nucFreqMap;
}
