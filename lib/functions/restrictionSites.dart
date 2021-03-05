import 'package:biokit/constants/strings.dart';
import 'package:biokit/functions/compStrand.dart';
import 'package:biokit/helpers/genCombos.dart';
import 'package:meta/meta.dart';

Map<String, List<Map<String, int>>> restrictionSites({
  @required nucSeq,
  @required seqType,
  int minSiteLen = 4,
  int maxSiteLen = 12,
  bool pos = false,
}) {
  String revCompSeq =
      compStrand(nucSeq: nucSeq, seqType: seqType, reversed: true);
  List<String> origSeqCombos = genCombos(seq: nucSeq);
  Iterable<String> seqIter = origSeqCombos
      .where((seq) => (seq.length >= minSiteLen) && (seq.length <= maxSiteLen));
  List<String> restSeqs = seqIter.toSet().toList();

  Map<String, List<Map<String, int>>> restSiteSeqs = {};
  List<Map<String, int>> restSiteLocations = [];

  for (String restSeq in restSeqs) {
    restSiteLocations = [];
    RegExp regRestExp = RegExp(restSeq);
    Iterable<RegExpMatch> matches = regRestExp.allMatches(revCompSeq);
    for (RegExpMatch match in matches) {
      //
      int startIdx = nucSeq.length - (restSeq.length + match.start);
      int endIdx = startIdx + restSeq.length;
      if (nucSeq.substring(startIdx, endIdx) == restSeq) {
        pos
            ? restSiteLocations
                .add({kStartPos: startIdx + 1, kEndPos: endIdx + 1})
            : restSiteLocations.add({kStartIndex: startIdx, kEndIndex: endIdx});
        restSiteSeqs[restSeq] = restSiteLocations;
      }
      ;
    }
  }
  return restSiteSeqs;
}
