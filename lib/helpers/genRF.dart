import 'package:biokit/constants/strings.dart';
import 'package:biokit/functions/compStrand.dart';
import 'package:biokit/helpers/nucSeqToAA.dart';
import 'package:biokit/helpers/reverseNucSeq.dart';
import 'package:meta/meta.dart';

/// Generates 6 ORFs, 3 forward, 3 reversed.
genRF({@required String nucSeq, @required String seqType}) {
  List<String> readingFrames = [];

  for (var i = 0; i < 3; i++) {
    readingFrames
        .add(nucSeqToAA(nucSeq: nucSeq, seqType: seqType, startIdx: i)[kAASeq]);

    readingFrames.add(nucSeqToAA(
        nucSeq: compStrand(nucSeq: nucSeq, seqType: seqType, reversed: true),
        seqType: seqType,
        startIdx: i)[kAASeq]);
  }
  return readingFrames;
}
