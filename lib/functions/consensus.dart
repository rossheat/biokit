import 'package:biokit/constants/lists.dart';
import 'package:biokit/constants/strings.dart';
import 'package:biokit/helpers/profileMat.dart';
import 'package:meta/meta.dart';

consensus(
    {@required List<String> nucSeqs,
    @required String seqType,
    bool rosPrint = false}) {
  Map<String, List<int>> profile =
      profileMat(nucSeqs: nucSeqs, seqType: seqType);

  String consensus = '';
  int seqLen = profile[profile.keys.first].length;

  for (var i = 0; i < seqLen; i++) {
    String mostFreqNuc = '';
    int mostFreqNucValue = 0;
    for (var nuc in seqType == kDNA ? dnaNucs : rnaNucs) {
      if (profile[nuc][i] > mostFreqNucValue) {
        mostFreqNuc = nuc;
        mostFreqNucValue = profile[nuc][i];
      }
    }
    consensus += mostFreqNuc;
  }

  return {'consensus': consensus, ...profile};
}
