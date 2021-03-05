import 'package:biokit/helpers/genRF.dart';
import 'package:biokit/helpers/rfToProtein.dart';
import 'package:meta/meta.dart';

List<String> nucSeqToProteins(
    {@required String nucSeq, @required String seqType, bool unique = false}) {
  List<String> rfs = genRF(nucSeq: nucSeq, seqType: seqType);

  List<String> proteins = [];
  for (var rf in rfs) {
    List<String> tempProteins = rfToProtein(aaSeq: rf);
    if (proteins != []) {
      proteins.addAll(tempProteins);
    }
  }

  proteins.sort((b, a) => a.length.compareTo(b.length));

  if (unique) {
    return proteins.toSet().toList();
  }
  return proteins;
}
