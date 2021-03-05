import 'package:biokit/constants/strings.dart';
import 'package:biokit/helpers/motifToRe.dart';
import 'package:meta/meta.dart';

Map<String, dynamic> matchMotif(
    {@required String seq,
    @required String motifSeq,
    overlap = true,
    bool pos = false}) {
  List<Map<String, dynamic>> matchData = [];
  Map<String, dynamic> matchMotifResultMap = {};
  String tempRegexMotif = proMotifToRe(motifSeq: motifSeq);

  RegExp regexMotif =
      overlap ? RegExp('(?=$tempRegexMotif)') : RegExp(tempRegexMotif);

  Iterable<RegExpMatch> allMatches = regexMotif.allMatches(seq);

  for (RegExpMatch match in allMatches) {
    pos
        ? matchData.add({
            kMatch: motifSeq,
            kStartPos: match.start + 1,
            kEndPos: match.start + motifSeq.length
          })
        : matchData.add({
            kMatch: motifSeq,
            kStartIndex: match.start,
            kEndIndex: match.start + motifSeq.length - 1
          });
  }

  matchMotifResultMap[kMatchCount] = allMatches.length;
  matchMotifResultMap[kMatchIndices] = matchData;
  return matchMotifResultMap;
}
