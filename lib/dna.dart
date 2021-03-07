import 'package:biokit/strings.dart';
import 'package:biokit/structs.dart';
import 'package:meta/meta.dart';
import 'package:biokit/nucleotides.dart';
import 'dart:math';

class DNA extends Nucleotides {
  DNA({required String seq}) : super(seq: seq, type: 'dna');

  /// Transcribes DNA to RNA.
  String transcribe() => this.seq.replaceAll('T', 'U');

  Map<String, List<Map<String, int>>> restrictionSites({
    int minSiteLen = 4,
    int maxSiteLen = 8,
  }) {
    String revComp = complementary(rev: true);

    List<String> origSeqCombos = combinations();
    Iterable<String> seqIter =
        origSeqCombos.where((seq) => (seq.length >= minSiteLen) && (seq.length <= maxSiteLen));
    List<String> restSeqs = seqIter.toSet().toList();

    Map<String, List<Map<String, int>>> restSiteSeqs = {};
    List<Map<String, int>> restSiteLocations = [];

    for (String restSeq in restSeqs) {
      restSiteLocations = [];
      RegExp regRestExp = RegExp(restSeq);
      Iterable<RegExpMatch> matches = regRestExp.allMatches(revComp);
      for (RegExpMatch match in matches) {
        int startIdx = this.seq.length - (restSeq.length + match.start);
        int endIdx = startIdx + restSeq.length;
        if (this.seq.substring(startIdx, endIdx) == restSeq) {
          restSiteLocations.add({kStartIndex: startIdx, kEndIndex: endIdx});
          restSiteSeqs[restSeq] = restSiteLocations;
        }
        ;
      }
    }
    return restSiteSeqs;
  }

  double tranRatio({required DNA oSeq}) {
    if (this.len != oSeq.len) {
      throw ('Unequal Sequence Lengths Error.');
    }

    int transitionCount = 0;
    int transversionCount = 0;

    this.seq.split('').asMap().forEach((idx, nuc) {
      if (nuc != oSeq.seq[idx]) {
        if (Structs.dnaTransitions.contains(nuc + oSeq.seq[idx])) {
          transitionCount++;
        } else if (Structs.dnaTransversions.contains(nuc + oSeq.seq[idx])) {
          transversionCount++;
        }
      }
    });

    return transitionCount / transversionCount;
  }

  static DNA random({required int len}) {
    Random _rand = Random();
    String dnaNucsStr = Structs.dnaNucs.join();
    String seq = String.fromCharCodes(
      Iterable.generate(
        len,
        (_) => dnaNucsStr.codeUnitAt(
          _rand.nextInt(dnaNucsStr.length),
        ),
      ),
    );
    return DNA(seq: seq);
  }
}
