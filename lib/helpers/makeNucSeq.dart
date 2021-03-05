import 'dart:math';
import 'package:meta/meta.dart';
import '../constants/lists.dart';
import '../constants/strings.dart';
import 'validateSeqType.dart';

String makeNucSeq({@required int seqLen, @required String seqType}) {
  String validSeqType = validateSeqType(seqType: seqType);

  Random _rand = Random();

  if (validSeqType == kDNA) {
    String dnaNucsStr = dnaNucs.join();
    return String.fromCharCodes(
      Iterable.generate(
        seqLen,
        (_) => dnaNucsStr.codeUnitAt(
          _rand.nextInt(dnaNucsStr.length),
        ),
      ),
    );
  }

  String rnaNucsStr = rnaNucs.join();
  return String.fromCharCodes(
    Iterable.generate(
      seqLen,
      (_) => rnaNucsStr.codeUnitAt(
        _rand.nextInt(rnaNucsStr.length),
      ),
    ),
  );
}
